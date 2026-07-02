# P2-6: FEM 力学ソルバー (線形弾性 + Maxwell 緩和)

## 目的

現行の `mechanics.cpp` は Maxwell 構成則の「セル毎更新式」のみで、
変位場もひずみも解いていない。四面体 P1 要素の線形弾性 FEM を実装し、
熱失配・真性膜応力を荷重として応力場を解けるようにする。
応力依存の拡散/酸化 (P3) の基盤。2 コミットに分割してよい
(コミット 1: K 組立 + patch test、コミット 2: 荷重・BC・緩和・API)。

## 現状コード

- `src/mechanics.cpp` — `ViscoElasticParams` (E, nu, tau; Si 既定
  130 GPa / 0.28)、`maxwell_update(sigma6, deps6, params, dt)`
  (Voigt 6 成分)。**流用する** (解いた弾性応力に緩和を適用)
- `src/sparse.cpp` — `cg_ilu0` (S-1 で並列化済み) を K u = f に使う
- `include/cprocess/mesh.hpp` — nodes/cells/cell_region/cell_vol、
  境界パッチ (`find_patch("zmin")` 等)
- P1-10 ParamDB (前提): 材料毎パラメータのキーに使う

## 実装手順

### コミット 1: 剛性行列と patch test

1. `include/cprocess/fem.hpp` + `src/fem.cpp` 新設:
   ```cpp
   struct FemProblem {
     CSR K;                        // 3*n_nodes 自由度、スカラー CSR
     std::vector<double> f;       // 荷重ベクトル
     std::vector<int> dof_fixed;  // 固定 dof (0/1)、サイズ 3*nn
     std::vector<double> dof_val; // 固定値 (通常 0)
   };
   ```
2. 定ひずみ四面体 (CST tet) の B 行列: 頂点座標から形状関数勾配
   `∇N_a = grad_a` を求める閉形式を仕様に明記する:
   - 体積 V = signed_vol(x0,x1,x2,x3)
   - `grad_0 = (x2−x1)×(x3−x1) / (6V)` (a=0 の対向面の面積ベクトル /3V、
     符号は内向き正になるよう順序に注意 — 4 頂点の巡回で
     `grad_a = (−1)^a (x_{b}−x_{d})×(x_{c}−x_{d}) / (6V)` 形式を採用し
     検算条件 Σ_a grad_a = 0 を組立時に assert)
   - B (6×12, Voigt): 行 (xx,yy,zz,yz,xz,xy)、列はノード a の (ux,uy,uz)。
     `B[0][3a]=gx, B[1][3a+1]=gy, B[2][3a+2]=gz, B[3][3a+1]=gz,
      B[3][3a+2]=gy, B[4][3a]=gz, B[4][3a+2]=gx, B[5][3a]=gy, B[5][3a+1]=gx`
   - D (6×6 等方): λ = Eν/((1+ν)(1−2ν)), µ = E/(2(1+ν));
     D の左上 3×3 = λ + 2µδ、右下対角 = µ
   - 要素剛性 k_e = V · Bᵀ D B (12×12)
3. アセンブリ: dof id = 3*node + comp。パターン構築はノード隣接
   (`MeshTopology::node_adj` + 自身) から 3×3 ブロックを展開して
   スカラー CSR に。組立は要素ループでスカラー加算
   (逐次で可 — 並列化は PA 系後続)
4. Dirichlet BC は **消去法ではなく対角固定**: 固定 dof の行を
   単位行列化 (対角 1、非対角 0、rhs = 固定値)、対応する列の寄与を
   rhs から差し引いてから列も 0 化 — SPD 保存の標準手順を仕様に明記
5. **Patch test** (テスト 1 で使用): 全境界ノードに線形変位場
   u = A·x (任意の 3×3 行列 A、仕様では A = [[1e-4,2e-5,0],[0,-3e-5,1e-5],[0,0,5e-5]]) を
   Dirichlet 指定 → 内部ノードの解も u = A·x を再現するはず

### コミット 2: 荷重・プロセス API・緩和

6. 荷重 (固有ひずみ法): 材料 m ごとに
   - 熱失配: ε0 = α_m·ΔT·I (体積ひずみ)、α は ParamDB キー
     `"mech.alpha.silicon"`=2.6e-6, `"mech.alpha.oxide"`=0.5e-6,
     `"mech.alpha.nitride"`=3.3e-6 [1/K]
   - 真性膜応力: σ0_m (ParamDB `"mech.sigma0.nitride"`=1e10 dyn/cm²
     ≈1 GPa 引張等) → ε0 = D⁻¹σ0
   - 要素荷重 f_e += V · Bᵀ D ε0
7. BC: zmin 全固定 (u=0)、側面は法線成分のみ固定 (roller)。
   パッチ判定は `find_patch("xmin")` 等 + 境界ノード集合
8. 解いた u から要素応力 σ = D(B u − ε0) を計算し、
   `maxwell_update` で dt_s 分の粘弾性緩和を適用 (材料の tau は
   ParamDB `"mech.tau.oxide"` 等; Si は tau=∞=緩和なし)。
   結果を `st.fields["sxx"].."sxy"` (6 成分) に格納
9. API 三層:
   - `proc::mechanics(SimState& st, double temp_k, double dt_s,
     std::ostream* log)` — ΔT = temp_k − 300K 基準
   - pybind `proc_mechanics(state, temp_k, dt_s)`
   - `Simulation.mechanics(temp, time)` (°C/min 変換)
   - デッキ `mechanics temp=1000C time=10min`

## テスト仕様

`tests/test_fem.cpp` 新設 (CMakeLists foreach 登録):

1. **Patch test**: 2×2×2 箱メッシュ、線形変位場 BC → 内部ノード変位の
   誤差 < 1e-10 (相対)、要素応力が一様 (全要素の σ_xx 差 < 1e-8 相対)
2. **一様熱膨張**: 単一材料 Si、roller 側面 + zmin 固定、ΔT=900K →
   偏差応力 (von Mises) < 1 MPa (自由膨張で応力ゼロ; 拘束方向の
   直応力は許容 — von Mises で判定する理由を仕様に記載)
3. **バイマテリアル**: Si 基板 + oxide 膜 (deposit で作成)、ΔT=900K →
   膜の面内応力 σ_xx < 0 (圧縮)、基板表面直下 σ_xx > 0 (引張)、
   膜応力の大きさが 10–1000 MPa のオーダー
4. **Maxwell 緩和**: テスト 3 の応力に tau=60s の oxide で dt=600s →
   膜応力が初期値の < 20% に減衰
5. 既存テスト全 PASS

Python: `Simulation.mechanics` スモーク (sxx フィールド生成、finite、
バイマテリアルの符号) を `python/test_comprehensive.py` に追加。

## 完了条件 (DoD)

- [ ] コミット 1: K 組立 + patch test PASS (`P2-6a` を含むメッセージ)
- [ ] コミット 2: 荷重/BC/緩和/API 三層 + テスト 2-5 PASS (`P2-6b`)
- [ ] 応力 6 成分が VTU 出力に載る (fields 経由で自動)
- [ ] 全既存テスト PASS

## やらないこと

- 大変形・幾何非線形、動的解析
- 応力依存の拡散/酸化への結線 (P3)
- 酸化の粘性流モデルとの統合 (P2-4 側で参照)
- 並列アセンブリ (後続 PA タスク)
