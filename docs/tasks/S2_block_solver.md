# S-2: ブロック連立ソルバー (ブロック CSR + ブロック Jacobi BiCGSTAB)

## 目的

P2-1 (フル点欠陥モデル) では dopant + I + V の複数種がセルごとに強く結合した
連立方程式になる。現行の CSR はスカラー行列のみで、種間結合は Picard の外側
反復に押し出すしかなく、硬い反応項では収束しない。本タスクでは各エントリが
密な nb×nb ブロック (nb = 連成種数、実行時サイズ) であるブロック CSR 構造と、
ブロック Jacobi 前処理付き BiCGSTAB を追加する。物理 (P2-1 本体) は含まない —
ソルバー基盤のみ。内部ライブラリ変更のため Python 表面 (pybind/Simulation) は
一切不要。

## 現状コード

- `include/cprocess/sparse.hpp` — `CSR { int n; ptr, col; val; }` (スカラー)、
  `Precond = std::function<void(const vector<double>&, vector<double>&)>`、
  ソルバー: `cg` / `bicgstab` (Precond 受け)、`*_jacobi` / `*_ilu0` ラッパー
- `src/sparse.cpp` — `bicgstab(A, b, x, rtol, maxit, psolve)`:
  約 60 行のループ。ベクトル演算はすべて素朴な `#pragma omp parallel for`
  ループ (`dotv`, `axpy`, 各種 elementwise) で、CSR に依存するのは
  `A.mul(...)` の 4 箇所のみ
- ブロック行列は存在しない

## 実装手順

1. `include/cprocess/sparse.hpp` に追加:
   ```cpp
   // Block CSR: square block matrix of n block-rows; each stored entry is a
   // dense nb x nb block (row-major). val.size() == ptr[n] * nb * nb.
   struct BCSR {
     int n = 0;    // block rows
     int nb = 0;   // block size (number of coupled species)
     std::vector<int> ptr;      // size n+1 (block entries)
     std::vector<int> col;      // size nnz_blocks, sorted within a row
     std::vector<double> val;   // nnz_blocks * nb*nb, block k at val[k*nb*nb]

     // y <- A x, with x/y of length n*nb (cell-major: x[i*nb + s]).
     void mul(const std::vector<double>& x, std::vector<double>& y) const;
     int find(int row, int c) const;  // block slot index, -1 if absent
   };

   // Builds a BCSR sharing the sparsity pattern of a scalar CSR (val zeroed).
   BCSR bcsr_from_pattern(const CSR& scalar_pattern, int nb);

   // Block-Jacobi preconditioned BiCGSTAB on a BCSR system.
   SolveResult bicgstab_bjacobi(const BCSR& A, const std::vector<double>& b,
                                std::vector<double>& x, double rtol, int maxit);
   ```
2. `src/sparse.cpp` に `BCSR::mul` を実装 (ブロック行並列):
   ```cpp
   void BCSR::mul(const std::vector<double>& x, std::vector<double>& y) const {
     y.assign(static_cast<std::size_t>(n) * nb, 0.0);
   #pragma omp parallel for schedule(static)
     for (int i = 0; i < n; ++i) {
       double* yi = &y[static_cast<std::size_t>(i) * nb];
       for (int k = ptr[i]; k < ptr[i + 1]; ++k) {
         const double* B = &val[static_cast<std::size_t>(k) * nb * nb];
         const double* xj = &x[static_cast<std::size_t>(col[k]) * nb];
         for (int r = 0; r < nb; ++r) {
           double s = 0;
           for (int c = 0; c < nb; ++c) s += B[r * nb + c] * xj[c];
           yi[r] += s;
         }
       }
     }
   }
   ```
   `BCSR::find` はスカラー `CSR::find` と同一の二分探索 (ブロックスロットを返す)。
3. `bcsr_from_pattern`: `ptr`/`col` をコピー、`n = scalar_pattern.n`、`nb` 設定、
   `val.assign(ptr[n] * nb * nb, 0.0)`。パターンが空 (`n==0`) なら空 BCSR を返す。
4. 密 nb×nb の部分ピボット付き in-place LU を無名 namespace に実装
   (nb は実行時値だが 3〜5 想定。固定サイズ特殊化は不要):
   ```cpp
   // In-place LU with partial pivoting of the nb x nb matrix `a` (row-major).
   // piv[k] records the row swapped into position k. Returns false if a pivot
   // is exactly zero (singular block).
   bool lu_factor_dense(double* a, int* piv, int nb);
   // Solves L U x = b using the factors/pivots; b is overwritten with x.
   void lu_solve_dense(const double* a, const int* piv, int nb, double* b);
   ```
   `lu_factor_dense` のアルゴリズム: 列 k = 0..nb-1 について
   (a) `|a[i*nb+k]|` 最大の行 i>=k を piv[k] とし行 k と swap、
   (b) `a[k*nb+k] == 0` なら false、
   (c) i > k について `a[i*nb+k] /= a[k*nb+k]`、
   `a[i*nb+j] -= a[i*nb+k]*a[k*nb+j]` (j > k)。
   `lu_solve_dense`: piv 順に b を並べ替えつつ前進代入 (単位下三角)、後退代入。
5. ブロック Jacobi 前処理: setup で各ブロック行 i の対角ブロック
   (`find(i, i)`) をコピーして `lu_factor_dense`。factor 失敗ブロックは
   単位行列扱い (apply が恒等)。apply は
   `z[i*nb..] = (D_ii)^{-1} r[i*nb..]` をブロック行並列で。
   `struct BlockJacobi { int n, nb; std::vector<double> lu; std::vector<int> piv;
   std::vector<char> ok; void setup(const BCSR&); void apply(const
   std::vector<double>& r, std::vector<double>& z) const; }` として
   sparse.cpp の無名 namespace でよい (ヘッダ公開不要)。
6. `bicgstab_bjacobi` の本体: **既存 `bicgstab(CSR, ...)` をテンプレート化
   しない**。既存関数の約 60 行 (r0, p, v, s, t, ph, sh のループ) をコピーし、
   `A.mul` を `BCSR::mul` に、`psolve` を `BlockJacobi::apply` に、ベクトル長を
   `n_total = A.n * A.nb` に置き換えた新関数として書く。`dotv` / `axpy` /
   elementwise ループは長さが変わるだけでそのまま流用できる (既に素朴な
   ループなので変更不要)。break 条件・収束判定 (`rn <= rtol * bnorm`)・
   `SolveResult` の埋め方は既存 bicgstab と一字一句同じ挙動にする。

## テスト仕様

新規 `tests/test_block_solver.cpp` (main + `test_util.hpp` の assert パターン)。
`CMakeLists.txt` の `foreach(t ...)` に `block_solver` を追加。

1. **nb=1 一致テスト**: 50×50 の 2D 5 点ラプラシアン (S-1 テストの
   `laplacian2d` 相当をテスト内で組み立て、n=2500) について、
   (a) スカラー `bicgstab_jacobi(A, b, x1, 1e-12, 2000)`、
   (b) `bcsr_from_pattern(A, 1)` に同じ値を詰めて
   `bicgstab_bjacobi(B, b, x2, 1e-12, 2000)`。
   両者 converged、かつ `max_i |x1[i] - x2[i]| < 1e-8`。
   (nb=1 のブロック Jacobi は点 Jacobi と同一なので数学的に同じ反復列)
2. **2 種連成 製造解テスト**: 8×8 グリッド (n=64 ブロック行、nb=2)。
   スカラーの 2D 5 点ラプラシアンパターンから `bcsr_from_pattern(P, 2)` を作り、
   対角ブロック = `[[4, -0.5], [-0.5, 4]]`、非対角ブロック = `[[-1, 0], [0, -1]]`
   を全エントリに設定。b は `b[i*2+s] = 1.0 + 0.01*(i*2+s)` の適当な
   非自明ベクトル。`bicgstab_bjacobi(B, b, x, 1e-12, 1000)` が収束し、
   テスト内で組んだ **密行列 (128×128) の Gauss 消去参照解** との
   `max |x - x_ref| < 1e-8`。密参照はテスト内の素朴な部分ピボット消去でよい。
3. **前処理の優位性**: テスト 2 と同じ系を
   (a) `bicgstab_bjacobi` (ブロック Jacobi)、
   (b) 同じ BCSR に対し点 Jacobi 前処理 (対角**スカラー** `1/A[i*nb+s, i*nb+s]`
   のみ; テスト内に ~15 行の点 Jacobi 版 apply を書いて `bicgstab_bjacobi` の
   コピーに差すか、BlockJacobi の対角外成分を 0 にした BCSR で代用) で解き、
   `iters_block < iters_point` を assert する。
   実装簡略のため「非対角成分を 0 にした対角ブロックのみで factor した
   BlockJacobi」を点 Jacobi の代用としてよい (数学的に同一)。
4. **スレッド決定性**: テスト 2 を `omp_set_num_threads(1)` と `(4)` で実行し、
   解の相対差 < 1e-9 (dot のリダクション順序が変わるため bit 一致は要求しない)。
   `#ifdef _OPENMP` ガードは test_solver.cpp / mc_implant.cpp のスタブパターン踏襲。

## 完了条件 (DoD)

- [ ] `BCSR` / `bcsr_from_pattern` / `bicgstab_bjacobi` が sparse.hpp/cpp に追加
- [ ] 既存のスカラー `bicgstab` は無変更 (diff が付かないこと)
- [ ] `tests/test_block_solver.cpp` の 4 テストが PASS、CMakeLists 登録済み
- [ ] 既存全テスト PASS (`ctest --test-dir build`)
- [ ] OpenMP 無効ビルドでもコンパイル・PASS
- [ ] コミットメッセージに `S-2` を含める

## やらないこと

- ブロック ILU(0) 前処理 (ブロック Jacobi で開始; 収束不足が実測されたら別タスク)
- 既存スカラー bicgstab のテンプレート化・共通化 (可読性優先で複製する)
- P2-1 の物理 (I/V 輸送方程式の組立) — 本タスクは行列構造とソルバーのみ
- CG 版 (連成行列は非対称なので BiCGSTAB のみ)
- DiffusionSolver への組込み・Python バインディング
