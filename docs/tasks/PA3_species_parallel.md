# PA-3: 種ループの並列化

## 目的

`DiffusionSolver::run()` の Picard 反復内では、共有量 (nni, dcell) を計算した
後の「種ごとの組立 + 線形ソルブ」が種間で独立している。多種連成
(B+P+As+Sb+I+V) では種数分の逐次ソルブが壁時間を支配するため、
種ループを OpenMP で並列化する。小さいメッシュ × 多種のケースで有効。

## 現状コード

- `src/diffusion.cpp` — `run()`: Picard ループ内の
  `for (int s = 0; s < ns; ++s) { assemble(...); cg_ilu0(A_, ...); }`。
  `A_.val`, `rhs`, `x`, `grad` は**メンバ/ローカルの共有バッファ**であり、
  そのままでは並列化するとデータ競合する
- `assemble()` は `A_.val` と引数 `rhs`/`grad` に書く。`A_` のパターン
  (ptr/col) と `diag_`/`fslot_`/`fg_`/`face_colors_` は読み取りのみ
- 内側の SpMV/AXPY/組立には既に OpenMP が入っている (S-1, 従来分)。
  **ネスト並列は既定で無効** (OpenMP default: max-active-levels=1) なので、
  外側を並列にすると内側は自動的に逐次実行になる

## 実装手順

1. `include/cprocess/diffusion.hpp`:
   - `DiffuseOpts` に `int species_parallel = 0;  // 0=auto, 1=on, -1=off` を追加
   - private に作業領域構造体を追加:
     ```cpp
     struct SolveWorkspace {
       std::vector<double> Aval, rhs, x;
       std::vector<Vec3> grad;
     };
     ```
2. `src/diffusion.cpp` — `run()` の変更:
   - Picard 反復の種ループ直前で有効判定:
     ```cpp
     const bool sp = (o.species_parallel > 0) ||
                     (o.species_parallel == 0 && ns >= 3 &&
                      nc < 50000);
     ```
     ヒューリスティクスの根拠をコメントで明記する
     (大メッシュでは内側 SpMV 並列の方が効率的; 小メッシュ多種では
     外側並列が勝つ)
   - `sp == false`: 現行コードパスを一切変更しない (バイト等価)
   - `sp == true`: ワークスペースプール `std::vector<SolveWorkspace> ws(ns)`
     を確保し (Aval は `A_.val.size()`、他は nc)、
     ```cpp
     #pragma omp parallel for schedule(dynamic, 1)
     for (int s = 0; s < ns; ++s) { ... }
     ```
     ループ本体では `assemble` を**ワークスペース版**で呼ぶ:
     `assemble_into(dcell[s], cold[s], bcface[s], c, dt, 0.0, o.nonortho,
                    ws[s].Aval, ws[s].rhs, ws[s].grad)` —
     既存 `assemble()` を `A_.val` の代わりに任意の val バッファへ書く
     形にリファクタし (`assemble()` は `assemble_into(..., A_.val, ...)` へ
     委譲)、ソルブは CSR のシャローコピー
     `CSR As{A_.n, A_.ptr, A_.col, ws[s].Aval}` (ptr/col はコピーで可、
     nc が小さい前提) に対して `cg_ilu0(As, ws[s].rhs, ws[s].x, ...)`
   - maxrel の集約は reduction(max:) または種ごとの部分値を配列に置き
     ループ後に max を取る (後者を指定: 決定的)
   - 例外: 並列リージョン内で throw しない — ソルバー不収束は
     種ごとの `failed[s]` フラグに記録し、ループ後に集約して throw
3. `run_ted()` は対象外 (I の解が種に先行依存するため)。コメントで明記

## テスト仕様

`tests/test_diffusion.cpp` に追加 (または `tests/test_species_parallel.cpp`
新設して CMakeLists の foreach に登録):

1. **等価性**: 8×8×8 箱メッシュ、B/P/As/Sb の 4 種を init + implant_gauss し、
   `species_parallel=-1` と `=1` で同一条件 diffuse (10min/1000C)。
   全種の全セルで相対差 < 1e-12 (各種のソルブは完全独立なので
   丸め順序も不変のはず — bit 一致が期待できるが判定は 1e-12 とする)
2. **auto 判定**: ns=1 のとき sp==false のパスを通る (species_parallel=0 で
   1 種 diffuse が現行テストと同一結果)
3. **失敗集約**: lin_maxit=1 で全種不収束 → throw されること
4. 既存の diffusion/ted/integration テストが全て PASS

Python: `python/test_comprehensive.py` に 4 種同時 diffuse の
スモークテスト (結果 finite、dose 保存) を追加。
DiffuseOpts.species_parallel の pybind 公開
(`.def_readwrite("species_parallel", ...)`) も行う。

## 完了条件 (DoD)

- [ ] `species_parallel` オプション実装、auto ヒューリスティクス動作
- [ ] off 時は現行コードパス不変 (回帰なし)
- [ ] 等価性テスト (on vs off < 1e-12) PASS
- [ ] 全既存テスト PASS
- [ ] コミットメッセージに `PA-3` を含める

## やらないこと

- run_ted の種並列化 (I 依存があるため)
- ネスト並列 (内側 OpenMP との同時有効化)
- ILU 分解の種間共有 (行列は種ごとに異なる)
- ns や nc によるヒューリスティクスの自動チューニング
