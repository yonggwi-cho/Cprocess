# S-4: AMG 前処理 (smoothed aggregation, 2 レベル)

## 目的

ILU(0) 前処理の CG 反復数はセル数とともに増大し、10⁶ セル級の 3D で律速に
なる。格子非依存の収束を狙い、smoothed-aggregation AMG を前処理として追加
する。まず **2 レベル V サイクル** を実装し (これが本タスクの合格範囲)、
再帰的マルチレベル化は明示のストレッチ項目とする。内部ライブラリ変更のため
Python 表面は不要。

## 現状コード

- `include/cprocess/sparse.hpp` — `CSR`、`Precond` (std::function)、
  `cg(A, b, x, rtol, maxit, psolve)` は任意の Precond を受ける → AMG は
  この Precond として差すだけで済む
- `src/sparse.cpp` — `cg_ilu0` / `bicgstab_ilu0` が「factor して Precond
  ラムダで包む」パターンの手本
- AMG は存在しない

## 実装手順

すべて `sparse.hpp` / `sparse.cpp` に追加 (新ファイル不要。肥大が気になる
場合は `src/amg.cpp` に分離し CMakeLists の `add_library` に追加してよい)。

1. **公開インターフェース** (`sparse.hpp`):
   ```cpp
   // Smoothed-aggregation AMG, two-level V-cycle. Usable as a Precond for
   // the existing cg/bicgstab cores.
   struct AMG {
     void setup(const CSR& A);   // build aggregates, P, R=P^T, Ac
     void apply(const std::vector<double>& r, std::vector<double>& z) const;
     int n_aggregates() const { return nagg_; }
    private:
     CSR A_;          // fine matrix (copy)
     CSR P_, R_;      // prolongation (n x nagg) and restriction (= P^T)
     CSR Ac_;         // coarse matrix R A P (nagg x nagg)
     std::vector<double> dinv_;       // fine 1/a_ii for the Jacobi smoother
     std::vector<double> coarse_lu_;  // dense LU of Ac when nagg < 200
     std::vector<int> coarse_piv_;
     bool coarse_dense_ = false;
     int nagg_ = 0;
   };
   SolveResult cg_amg(const CSR& A, const std::vector<double>& b,
                      std::vector<double>& x, double rtol, int maxit);
   ```
   注: `CSR` は正方前提のメンバ (`n`) しか持たないが、`mul` は
   `ptr/col/val` しか見ないので長方行列も `n = 行数` として扱える。
   P_ は n 行 (各行 1〜数 nnz)、R_ は nagg 行。`R_ = transpose(P_)` を作る
   ヘルパー `CSR csr_transpose(const CSR& A, int ncols)` を無名 namespace に
   実装 (カウント → prefix sum → scatter の標準 2 パス)。
2. **強結合グラフと集約** (setup 内):
   - 強結合判定: 非対角 (i,j) が `|a_ij| > θ * sqrt(a_ii * a_jj)`、θ = 0.08。
     a_ii ≤ 0 の行は自分だけの集約にする (堅牢化)
   - 貪欲集約: `agg[i] = -1` で初期化し、i = 0..n-1 順に
     (a) i が未割当かつ i の強結合近傍がすべて未割当なら、新しい集約 g を
     作り i と強結合近傍全員を g に入れる (root-node aggregation)。
     (b) 第 1 パス後に残った未割当ノードは、強結合近傍のうち既存集約に
     属すものがあればその (最初に見つかった) 集約へ吸収、なければ単独集約。
     決定的 (ノード番号順) にすること
   - `nagg_` = 集約数
3. **tentative prolongator と平滑化**:
   - P0: `P0[i, agg[i]] = 1.0` (piecewise constant; 各行 nnz=1 の CSR)
   - ρ̂ の推定: D^{-1}A に対する電力法 10 反復。初期ベクトルは
     `v[i] = 1.0 + 1e-3 * (i % 7)` の決定的なもの、各反復で
     `v <- D^{-1} (A v)` して `ρ̂ = ||v_new|| / ||v_old||` を最後の反復で採用
   - `ω = 4.0 / (3.0 * ρ̂)`、`P = (I - ω D^{-1} A) P0`。
     P0 の行 i の nnz は 1 (列 agg[i]) なので、P の行 i の列集合は
     `{agg[j] : a_ij ≠ 0}`。実装: 行 i ごとに小さな (col → val) の
     線形探索付き小配列 (行 nnz は最大でも A の行 nnz) で
     `P[i, agg[j]] += (δ_ij - ω * dinv[i] * a_ij) * 1.0` を蓄積し、
     列ソートして CSR 化
4. **粗行列 Ac = R A P** (疎三重積):
   - まず `AP = A * P` を行ごとに計算: 行 i について
     `std::unordered_map<int, double> acc;` に
     `for k in A.row(i): for l in P.row(A.col[k]): acc[P.col[l]] += A.val[k]*P.val[l]`
     を蓄積 → キーをソートして CSR 行に書き出す (行並列 OpenMP 可、
     行ごとに独立バッファ)
   - `Ac = R * AP` を同じハッシュ蓄積で計算 (R の行 = 集約、行数 nagg_)。
     Ac の行内列は必ずソートすること (`CSR::find` の二分探索前提)
5. **スムーザと粗解法**:
   - pre/post スムーザ: 減衰 Jacobi 1 スイープ、ω_J = 2/3:
     `x <- x + ω_J D^{-1} (r - A x)` (apply 内では初期 x=0 なので
     pre は `x = ω_J D^{-1} r` に簡約できる)
   - 粗解法: `nagg_ < 200` なら setup 時に Ac を密化して部分ピボット LU を
     保存 (S-2 の `lu_factor_dense`/`lu_solve_dense` と同一アルゴリズム。
     S-2 未実装ならローカルに同じものを書く)、apply では代入のみ。
     `nagg_ >= 200` なら減衰 Jacobi 10 スイープで近似解
6. **apply (V(1,1) サイクル)**:
   ```
   x  = ωJ D^{-1} r                    // pre-smooth (x0 = 0)
   rf = r - A x                        // fine residual
   rc = R rf                           // restrict
   ec = coarse_solve(Ac, rc)           // dense LU or Jacobi sweeps
   x += P ec                           // prolongate + correct
   x += ωJ D^{-1} (r - A x)            // post-smooth
   z  = x
   ```
   すべて既存の `CSR::mul` + elementwise ループで書ける (OpenMP 済み)。
7. `cg_amg`: `AMG amg; amg.setup(A);` + Precond ラムダで `cg(...)` に委譲
   (`cg_ilu0` と同型)。
8. **ストレッチ (任意; DoD 外)**: `setup` を再帰化し、`nagg_ >= 200` の
   場合に粗レベルへさらに AMG を掛ける多レベル版。実装する場合は
   `std::unique_ptr<AMG> next_;` を持たせ、粗解法を `next_->apply` に差し
   替えるだけの構造にする。テスト基準は 2 レベル版と同じでよい。

## テスト仕様

新規 `tests/test_amg.cpp`、CMakeLists の foreach に `amg` を追加。
テスト内ヘルパー `CSR laplacian3d(int n)` (n³ の 7 点ステンシル、
対角 6、隣接 -1、Dirichlet 境界) を実装。

1. **反復数**: n=32 (32768 行) の 3D ラプラシアン、b は全要素 1。
   `cg_amg(A, b, x, 1e-10, 500)` が収束し、
   (a) **iters < 40**、(b) 同条件の `cg_ilu0` の反復数の **1/2 未満**。
2. **解の一致**: 上記の `cg_amg` 解と `cg_ilu0` 解の
   `max_i |x_amg[i] - x_ilu[i]| < 1e-8`。
3. **時間 (情報のみ)**: `std::chrono` で setup と solve の wall time を
   `printf` する。assert しない。
4. **集約数**: `n/20 <= amg.n_aggregates() <= n/3` (n = 32768)。
   7 点ステンシルの近傍集約なら典型 ~n/7〜n/15 に入る。
5. **スレッド決定性**: n=16 で 1 スレッドと 4 スレッドの `cg_amg` 解の
   相対差 < 1e-9 (集約・P 構築は決定的、dot リダクションのみ順序が変わる)。

## 完了条件 (DoD)

- [ ] `AMG::setup/apply` と `cg_amg` が追加され既存 Precond 機構で動く
- [ ] `tests/test_amg.cpp` の 5 項目 PASS (反復数 <40 かつ ILU0 の半分未満)
- [ ] 既存全テスト PASS、OpenMP 無効ビルドでも PASS
- [ ] setup が決定的 (同一入力 → 同一集約・同一 P) であること
- [ ] コミットメッセージに `S-4` を含める

## やらないこと

- DiffusionSolver への組込み (前処理の選択ロジック変更は別タスク。
  当面は ILU0 のままで、大規模時に切替える)
- W サイクル、K サイクル、Chebyshev/Gauss-Seidel スムーザ
- 非対称系向けの smoothed aggregation 調整 (CG/SPD 前提で開始)
- マルチレベル再帰 (ストレッチ扱い; 2 レベルで DoD を満たす)
- Python バインディング
