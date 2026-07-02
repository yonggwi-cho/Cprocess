# S-3: Newton-Krylov (JFNK) オプション

## 目的

`DiffusionSolver::run` の非線形性 (濃度依存拡散係数 D(c)、電場増速) は現在
Picard 反復で処理しているが、高濃度 (強非線形) では `max_picard` の予算を
使い切っても収束が遅い。Jacobian-free Newton-Krylov (JFNK) を代替オプション
として追加する: 残差 F(c) を既存 assemble 機構で評価し、Jacobian ベクトル積を
有限差分方向微分で近似、線形化系を GMRES で解く。P2-1/P2-2 の硬い反応項への
布石。内部ソルバー変更のため Python 表面は不要 (`DiffuseOpts` は proc:: 層から
渡されるが新フィールドはデフォルト値で後方互換)。

## 現状コード

- `src/diffusion.cpp` — `DiffusionSolver::run`: 時間ステップごとに Picard
  ループ (`for picard = 1..o.max_picard`)。各反復で
  (a) 全種の現在反復値から `nni` (n/ni) を計算、
  (b) 種 s ごとに `dcell[s]` (D(c)) を評価、
  (c) `assemble(dcell[s], cold[s], bcface[s], c, dt, 0.0, o.nonortho, rhs, grad)`
  で `A_` (メンバ) と `rhs` を組み、`cg_ilu0` で解く。
  収束判定は反復間の相対変化 `maxrel < o.picard_tol`
- `include/cprocess/diffusion.hpp` — `DiffuseOpts` (temp, time, dt,
  max_picard, picard_tol, lin_rtol, lin_maxit, ...)
- `src/sparse.cpp` — `gmres_jacobi(A, b, x, rtol, maxit, restart)`:
  約 110 行の restart 付き GMRES。`A.mul(t, r)` (2 箇所) と
  `inv_diag(A)` による右 Jacobi 前処理以外は行列に依存しない
- `ILU0` — `factor(A)` + `apply(x, y)` (前処理として使用可能)

## 実装手順

1. **GMRES の演算子化リファクタ** (`sparse.hpp` / `sparse.cpp`):
   ```cpp
   // Operator form: aop(x, y) computes y <- A x. Right-preconditioned with
   // psolve (y <- M^{-1} x); pass an identity lambda for none.
   using LinOp = std::function<void(const std::vector<double>&,
                                    std::vector<double>&)>;
   SolveResult gmres_op(const LinOp& aop, int n, const std::vector<double>& b,
                        std::vector<double>& x, double rtol, int maxit,
                        int restart, const Precond& psolve);
   ```
   実装: 既存 `gmres_jacobi` の本体を移し、
   (a) `A.mul(x, t)` → `aop(x, t)`、
   (b) `M[i] * v[i]` (Jacobi 適用、2 箇所: Arnoldi 内と解更新) →
   `psolve(v, tmp)` に置換、(c) `A.n` → 引数 `n`。
   `gmres_jacobi` は薄いラッパーに書き換える:
   `aop = [&](x,y){ A.mul(x,y); }`, `psolve = Jacobi(inv_diag(A))`。
   既存テストの gmres 結果が変わらないこと (数学的に同一の演算列)。
2. **DiffuseOpts 拡張** (`diffusion.hpp`):
   ```cpp
   bool use_newton = false;    // per-step nonlinear solve: JFNK instead of Picard
   double newton_rtol = 1e-8;  // ||F|| < newton_rtol * ||F0|| convergence
   ```
3. **残差評価ヘルパー** (`DiffusionSolver` の private メソッド):
   ```cpp
   // F(c) for one backward-Euler step of species s:
   //   F = A(c) * c - rhs(c)
   // where assemble(dcell(c), cold, bcface, c, dt, ...) builds A_ and rhs.
   // dcell(c) is re-evaluated from the trial c (nni fixed from the outer
   // Newton iterate's all-species state — see step 5).
   void residual(const std::vector<double>& dcell_c,
                 const std::vector<double>& cold,
                 const std::vector<double>& bcface,
                 const std::vector<double>& c, double dt, bool nonortho,
                 std::vector<double>& F);
   ```
   実装: `assemble(dcell_c, cold, bcface, c, dt, 0.0, nonortho, rhs_, grad_)`
   を呼んで `A_`/`rhs_` を作り、`A_.mul(c, F)` の後
   `F[i] -= rhs_[i]` (OpenMP for)。rhs_/grad_ は使い回しのメンバまたは
   呼び出し側スクラッチ。**注意**: assemble は A_ を上書きするため、
   residual 呼出し後の A_ は「その c で組んだ行列」になる — これを
   前処理凍結 (step 5) に利用する。
4. **JFNK の Jacobian-vector 積** (無名 namespace のヘルパーかラムダ):
   ```
   Jv ≈ (F(c + eps*v) - F(c)) / eps
   eps = sqrt(DBL_EPSILON) * (1 + ||c||_2) / ||v||_2    (||v||=0 なら Jv=0)
   ```
   F(c + eps*v) の評価では dcell も摂動後の c から再評価する
   (種 s の D は自種 c と nni に依存; nni は Newton 外側で凍結し
   自種依存分のみ厳密に微分する — Picard と同じ種分離近似)。
5. **Newton ループ** (`run` の Picard ループを置換する分岐):
   `o.use_newton == true` のとき、種 s ごとに:
   - (a) 外側 (種間) 結合は現行どおり: nni を全種の現在値から計算し凍結
     (Picard の外側 1 パスに相当。種ループの外側の反復回数・構造は
     現行の picard ループ枠をそのまま使い、内側の「1 回の線形ソルブ」を
     「1 回の Newton ソルブ」に置き換える)
   - (b) `c0 = c` (現在値) で dcell を評価、`residual(...)` で `F0` を計算、
     `normF0 = ||F0||_2`。`normF0 == 0` なら収束済み
   - (c) この時点の `A_` (= A(c0)) から `ILU0 ilu; ilu.factor(A_);` —
     Newton 反復中は**凍結** (再分解しない)
   - (d) Newton 反復 k = 1..10:
     * `gmres_op(Jop, nc, Fneg, delta, o.lin_rtol, o.lin_maxit, 30, iluP)` で
       J δ = -F を解く (`Fneg[i] = -F[i]`、`iluP` = ilu.apply、
       `Jop` = step 4 の有限差分積。Jop 内の F 評価 1 回 = assemble 1 回 +
       SpMV 1 回であることをコメントで明記)
     * `c += delta`、負値は 0 にクランプ (現行と同じ後処理)
     * dcell を新しい c で再評価し `residual` で F 更新、
       `||F|| < o.newton_rtol * normF0` で収束 break
   - (e) k = 10 まで到達して未収束なら
     `throw std::runtime_error("diffusion: Newton failed to converge in 10 "
     "iterations (||F||/||F0||=" + std::to_string(...) + ")")`
   - ログ (`verbosity>=1`) の `picard=` フィールドに Newton 反復数を出す
     (フォーマット互換のため同じ桁位置でよい; `newton=%d` に変えてもよいが
     既存文言に依存するテストがないことを確認して選ぶ)
6. `o.use_newton == false` (デフォルト) の経路は**一切変更しない**。

## テスト仕様

新規 `tests/test_newton.cpp`、`CMakeLists.txt` の foreach に `newton` を追加。
メッシュは `make_box_mesh` の 12×12×24 程度 (test_diffusion.cpp の使い方踏襲)。

1. **Picard との一致**: B の Gaussian 初期プロファイル
   (peak 1e19 cm^-3、既存 test_diffusion のセットアップ流用)、
   1000 °C 相当 T=1273.15 K、time=600 s。
   `use_newton=false` と `use_newton=true` で別々に走らせ、
   peak と dose (Σ c·V) の相対差がともに **< 1e-6**。
2. **高濃度での優位性**: peak 5e20 cm^-3 (強い n/ni 非線形 + 電場増速) の
   B アニール、`max_picard=8, picard_tol=1e-10` (Picard に厳しい設定)。
   ステップあたりの反復数を比較できるよう、`DiffuseOpts` に
   `mutable int* nl_iters = nullptr;`(実行後: 全ステップの非線形反復数合計)
   を追加して両経路で記録する。assert:
   **newton_iters < picard_iters** かつ newton の 1 ステップ平均 ≤ 5。
3. **停滞時の例外**: `use_newton=true, newton_rtol=1e-30` (到達不能な要求) で
   run が `std::runtime_error` を投げること (`try/catch` で捕捉し
   メッセージに "Newton" を含むこと)。
4. **gmres_op リグレッション**: `gmres_jacobi` の既存テスト
   (tests/test_solver.cpp) がラッパー化後もそのまま PASS すること。

## 完了条件 (DoD)

- [ ] `gmres_op` 追加、`gmres_jacobi` はラッパー化され既存テスト PASS
- [ ] `DiffuseOpts::use_newton / newton_rtol / nl_iters` 追加 (デフォルトで挙動不変)
- [ ] `use_newton=false` の結果が本タスク前と bit 一致 (コード経路無変更)
- [ ] `tests/test_newton.cpp` の 4 項目 PASS、CMakeLists 登録
- [ ] 既存全テスト PASS
- [ ] コミットメッセージに `S-3` を含める

## やらないこと

- `run_ted` への Newton 適用 (スコープ外; run のみ。TED は P2-1 置換時に統合)
- 解析的 Jacobian の組立 (JFNK の趣旨に反する)
- ラインサーチ / トラストリージョン等のグローバル化 (max 10 反復 + throw で十分)
- 種間 (nni) 結合の Newton 内取込み (S-2 ブロックソルバー後の課題)
- Python バインディング / proc:: API の変更
