# PA-2: NUMA/帯域最適化 (RCM 並べ替えによる行列バンド幅縮小)

## 目的

SpMV はメモリ帯域律速であり、セル番号付け (箱型 Kuhn 分割や Gmsh 由来) が
空間的に散らばっていると `x[col[k]]` のアクセスがキャッシュに乗らない。
Reverse Cuthill-McKee (RCM) で行列バンド幅を縮小し、拡散ソルバー内の
線形ソルブのキャッシュ局所性を改善する。行列レベルの並べ替えに限定し、
メッシュのセル番号自体は変えない。内部変更のため Python 表面は不要。

## 現状コード

- `src/diffusion.cpp` — `DiffusionSolver::build()`: 面走査から CSR パターン
  `A_` を構築し、`diag_[i] = A_.find(i,i)`、
  `fslot_[fi] = {A_.find(owner,neigh), A_.find(neigh,owner)}` を保存。
  `assemble()` は `A_.val[diag_[i]]` / `A_.val[fslot_[fi][0/1]]` に直接
  scatter する (面彩色並列)。`run`/`run_ted` は `cg_ilu0(A_, rhs, x, ...)`
- `include/cprocess/sparse.hpp` — `CSR`。並べ替えユーティリティは無い
- セル順はメッシュ由来 (mesh.cpp / box_mesh.cpp) で変更しない

## 実装手順

### 1. RCM ユーティリティ (sparse.hpp / sparse.cpp)

```cpp
// Reverse Cuthill-McKee ordering of the (symmetric-pattern) matrix A.
// Returns perm with perm[new] = old.
std::vector<int> rcm_order(const CSR& A);

// Symmetric permutation: B(new_i, new_j) = A(perm[new_i], perm[new_j]).
// Row columns of the result are sorted.
CSR permute(const CSR& A, const std::vector<int>& perm);
```

- `rcm_order` の手順:
  1. 擬似周辺ノード探索: ノード 0 から BFS し最遠ノード u を得る。
     u から再度 BFS し最遠ノード v を得る。v を開始ノードとする
     (レベル数が減らなければ v で確定; 2 回で打ち切る)
  2. Cuthill-McKee: v から BFS。各ノードの未訪問隣接を**次数昇順・同次数は
     ノード番号昇順**で queue に積む (決定的)。非連結成分が残る場合は
     未訪問の最小番号ノードから同じ手順を繰り返す
  3. 得られた順序を反転して返す (reverse)
- `permute`: `iperm[old] = new` を作り、行 new_i = 行 perm[new_i] の
  各列 c を `iperm[c]` に写像、`(col, val)` ペアを列でソートして CSR 化。

### 2. DiffusionSolver への組込み (組立時リマップ方式)

方針: **ソルブ前後で permute するのではなく、A_ 自体を最初から並べ替え
済みパターンで持つ**。assemble の scatter 先スロット (`diag_`, `fslot_`)
を build 時に並べ替え済み行列のスロットへリマップしておけば、assemble は
無変更のまま並べ替え済み行列を直接組む。ソルブは permuted 空間で行い、
rhs/x のみ出入りで並べ替える。

1. `DiffusionSolver` に private メンバを追加:
   ```cpp
   std::vector<int> perm_;   // perm_[new] = old (cell index)
   std::vector<int> iperm_;  // iperm_[old] = new
   std::vector<double> rhs_p_, x_p_;  // permuted scratch
   ```
2. `build()` の末尾 (現行パターン確定後) に:
   - 現行手順で組んだ**スカラーパターン** A0 (現行 A_) から
     `perm_ = rcm_order(A0)`、`iperm_` を計算
   - `A_ = permute(A0, perm_)` (val はゼロ)
   - `diag_[i]` を再計算: `diag_[old_i]` は「old セル i の対角が入る
     A_ (permuted) 内のスロット」= `A_.find(iperm_[i], iperm_[i])`。
     つまり **diag_ / fslot_ の添字規約は従来どおり old セル/面番号のまま**、
     指す先だけ permuted 行列のスロットにする:
     `diag_[i] = A_.find(iperm_[i], iperm_[i]);`
     `fslot_[fi] = {A_.find(iperm_[owner], iperm_[neigh]),
                    A_.find(iperm_[neigh], iperm_[owner])};`
   - これで `assemble()` は 1 行も変えずに permuted 行列へ scatter する。
     ただし assemble 内の `rhs[i]` は old 番号のままなので、**rhs も
     permuted で持つ必要はない** — ソルブ直前に写像する (次項)
3. `run` / `run_ted` の各線形ソルブ箇所 (`cg_ilu0(A_, rhs, x, ...)` と
   bicgstab フォールバック、計 4 箇所) をヘルパーに集約する:
   ```cpp
   SolveResult solve_permuted(const std::vector<double>& rhs,
                              std::vector<double>& x, double rtol, int maxit);
   ```
   実装: `rhs_p_[new] = rhs[perm_[new]]`, `x_p_[new] = x[perm_[new]]`
   (warm start) → `cg_ilu0(A_, rhs_p_, x_p_, ...)` (+bicgstab fallback、
   現行の throw 文言維持) → `x[perm_[new]] = x_p_[new]`。
4. **バンド幅の計測とログ**: build() で
   `bw = max_i max_{k in row i} |i - col[k]|` を A0 (並べ替え前) と A_
   (後) について計算し、`log_` があれば
   `[diffuse] RCM: bandwidth 12345 -> 678 (n=NNN)` を出力する。
5. `mask_` で全セル無効などの縮退 (n==0) では perm_ を恒等にして通す。

## テスト仕様

新規 `tests/test_rcm.cpp` + `tests/test_diffusion.cpp` への追加。
CMakeLists の foreach に `rcm` を追加。

1. **rcm_order/permute 単体**: 16³ box メッシュ相当の 3D 7 点ラプラシアン
   (n=4096) で、
   - perm が置換であること (ソートして 0..n-1)
   - `permute` 後の行列で `A x` を計算し、並べ替え前の
     `A0 x0` (x0[old]=x[iperm 対応]) と全要素一致 (bit レベル、
     行内の加算順序は列ソートで固定される — 一致しない場合 < 1e-14 相対)
   - **バンド幅が 1/2 以下**: `bw_after * 2 <= bw_before`
     (7 点ラプラシアンの自然順は bw = n_x*n_y = 256、RCM で大幅縮小)
2. **結果不変**: 既存 test_diffusion の B アニールケースを本タスク前後で
   比較する代わりに、テスト内で「RCM を恒等 perm に強制」できる必要はなく、
   物理的な検証で足りる: dose 保存 (相対 < 1e-9)・peak が既存テストの
   合格基準のまま PASS すること。加えて解の再現性として、同一入力で
   run を 2 回実行し全セル bit 一致すること。
   参照値との一致は **相対 < 1e-9** を基準にする (同じ演算の順序変更に
   よる丸め差のみ許容; ILU0 は行内順序が変わるため bit 一致は要求しない)。
   実装方法: 本タスク着手前に既存ケースの peak/dose を出力して定数として
   テストに焼き込み、適用後に相対 < 1e-9 で一致を assert する。
3. **反復数の非劣化**: 上記ケースの `SolveResult::iters` 合計 (run 全体の
   lin_iters、ログから or opts に計測ポインタを足さず戻り値を printf) が
   適用前の ±10% 以内 (RCM は ILU0 の品質を通常改善する)。判定は
   テスト 2 と同様、適用前の実測値を定数化して assert する。
4. **既存全テスト PASS** (`ctest`)。

## 完了条件 (DoD)

- [ ] `rcm_order` / `permute` が追加され単体テスト PASS
- [ ] DiffusionSolver が permuted 行列で組立・求解し、結果が従来比
      相対 < 1e-9 で一致
- [ ] 16³ でバンド幅 ≥2 倍縮小、ログに before/after が出る
- [ ] 反復数 ±10% 以内、既存全テスト PASS、CMakeLists 登録
- [ ] コミットメッセージに `PA-2` を含める

## やらないこと

- **メッシュセル番号の Hilbert/Morton 並べ替え**: セル index は
  fields / cell_region / cell_locator / VTK 出力 / MC の hits 等
  全コンシューマに波及するため、行列レベルの局所化で先に効果を取る。
  効果不足が実測されたら別仕様で
- NUMA first-touch 初期化・明示的スレッドピニング (計測してから)
- ILU0 のフィルイン順最適化 (RCM の副次効果で足りる想定)
- run_ted 以外のソルバー利用箇所 (mechanics 等) への適用
- Python バインディング
