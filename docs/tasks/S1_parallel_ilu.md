# S-1: ILU(0) 分解フェーズの並列化

## 目的

`ILU0::factor()` (src/sparse.cpp) は逐次 IKJ ループであり、拡散ソルバーの
並列スケーリングを制限している (実測 1.95×/4 コア。三角ソルブは既に
レベルスケジューリングで並列)。分解にも同じレベルスケジューリングを適用し、
逐次ボトルネックを解消する。

## 現状コード

- `src/sparse.cpp` — `ILU0::factor(const CSR& A)`:
  1. `lu = A` をコピー、`diag` を計算
  2. **逐次** IKJ ループで in-place 分解 (行 i は行 k < i (a_ik≠0) に依存)
  3. 分解後にレベル集合 `lvl_lo` / `lvl_up` を構築 (これは解のためのもの)
- `ILU0::apply()` — レベルごとに `#pragma omp parallel for` で前進/後退代入 (変更不要)
- 呼び出し元: `cg_ilu0` / `bicgstab_ilu0` (毎回 factor を呼ぶ)

**鍵となる事実**: IKJ 分解における行 i の依存先は「行 k < i かつ A の
パターンで a_ik ≠ 0 の行 k」。これは前進代入の依存グラフと同一なので、
`lvl_lo` と同じレベル構造を**分解にも**使える。

## 実装手順

1. `ILU0::factor()` の処理順を変更する:
   - (a) `lu = A`; `diag` 計算 (現状どおり)
   - (b) **先に** `lvl_lo` / `lvl_up` を構築する。レベルはパターンのみに
     依存するので分解前に計算できる (現行コードの level 構築ブロックを
     IKJ ループの前に移動するだけ)
   - (c) IKJ ループを `lvl_lo` のレベル順に置き換える:
     ```cpp
     for (const auto& level : lvl_lo) {
       const int m = static_cast<int>(level.size());
     #pragma omp parallel for schedule(dynamic, 16)
       for (int t = 0; t < m; ++t) {
         const int i = level[t];
         // 現行の「行 i の分解」本体をそのまま移す:
         // for (kk = lu.ptr[i]; ...; ++kk) { k = lu.col[kk]; if (k >= i) break; ... }
       }
     }
     ```
   - 同一レベル内の行は互いに依存しないため、行 i の更新が読む
     `lu.val[dk]`, `lu.val[pk]` (行 k < i) は前レベルで確定済み。
     行 i 自身の書込み (`lu.val[kk]`, `lu.val[jj]`) は行 i の範囲のみ → 競合なし
2. `schedule(dynamic, 16)` を使う (レベル内の行の作業量は nnz に依存して不均一)
3. 結果は**逐次版と bit 一致**する (浮動小数の演算順序が行内で変わらないため)。
   これをテストで保証する

## テスト仕様

`tests/test_solver.cpp` に追加 (既存 main の ILU0 テスト群の後):

1. **bit 一致テスト**: 3D ラプラシアン様の疎行列 (下記) を作り、
   `omp_set_num_threads(1)` で factor した `lu.val` と
   `omp_set_num_threads(4)` で factor した `lu.val` が
   **全要素 `==` (bit 一致)** すること
   - 行列: `laplacian1d(n)` では帯幅が狭くレベル並列が効かないため、
     50×50 の 2D 5 点ラプラシアン (n=2500, 行あたり最大 5 nnz) を
     テスト内で組み立てる (新ヘルパー `laplacian2d(int nx, int ny)`)
   - OpenMP 無効ビルドでは自動的に両者同一になるので、テストはそのまま通る
2. **収束性テスト**: 同じ 2D ラプラシアンで `cg_ilu0` が
   `rtol=1e-10, maxit=500` 以内に収束し、製造解との最大誤差 < 1e-6
3. **既存テスト**: `cg_ilu0` の三重対角 1 反復収束テスト等がそのまま PASS すること

`#ifdef _OPENMP` で `<omp.h>` を include し、無効時は
`omp_set_num_threads` をスキップ (mc_implant.cpp のスタブパターンを踏襲)。

## 完了条件 (DoD)

- [ ] `ILU0::factor` がレベル並列化され、逐次フォールバック (OpenMP 無効) でも動く
- [ ] bit 一致テスト・収束性テストを含む全 16 テストが PASS
- [ ] 拡散ベンチ (`python3 /tmp/perf.py` 相当、166k セル) で 4 コアの
      wall-clock が現状 (~15s) から**悪化していない** (改善の確認は任意)
- [ ] コミットメッセージに `S-1` を含める

## やらないこと

- Chow-Patel 反復 ILU (レベル法で十分。将来 GPU 化時に再検討)
- ILU 分解の再利用 (Picard 反復間のキャッシュ) — 別タスク
- AMG やその他前処理の追加
