# P2-4: 2D/3D 酸化 (酸化剤拡散 + 界面法線速度場 — bird's beak)

## 目的

P1-6 の oxidize は 1D ブランケット (extend + 再タグ) であり、窒化膜マスク
下の横方向侵入 (bird's beak) を表現できない。本タスクで、窒化膜マスクが
存在する場合に:

1. 酸化膜領域内の**酸化剤 (O₂/H₂O) 定常拡散-反応**を解き、
2. Si/SiO₂ 界面の各面での**局所成長速度** v_n を求め、
3. **ノード毎法線速度場**版の ALE で界面を動かし、酸化膜の体積膨張
   (Si 消費の 2.27 倍、外面は 1.27 倍分持ち上げ) をノード単位で実現する。

応力 (粘性流・応力遅延酸化) は**明示的にスコープ外** (P3 (e) 参照)。
**P1-6、M-4 (界面適合)、M-5 (repair_quality) の完了が前提。**
大規模タスクのため **2 コミットに分割してよい** (下記 Stage A/B)。

## 現状コード

- `src/process.cpp` — `proc::oxidize` (P1-6/P2-3 後: サブステップ化済みの
  1D 経路)。窒化膜は `region_material` の `"nitride"`。
- `src/diffusion.cpp` — `DiffusionSolver`: CSR 組立 (`build()`/`assemble()`)、
  境界面透過率 `gb = |S|²/(S·(cf−cP))`、`cg_ilu0`。**セル mask で解領域を
  制限できる** (mask=1 のセルのみ解く) — 酸化膜サブメッシュ相当は
  「mask = oxide セルのみ 1」で実現できるが、Robin BC (界面反応) が
  ないため専用の小ソルバーを書く (下記)。P1-4 完了後は偏析の
  境界面透過率実装が Robin 型 BC の手本になる。
- `src/ale_mover.cpp` — `ale_move(m, topo, node_normals, v_n /*scalar*/, dt,
  mask_fn, smooth_iters)`: **スカラー速度のみ**。反転ガード +
  内部 Laplacian 平滑化。`compute_node_normals` は境界/界面ノードの
  面積重み法線を返す。
- `include/cprocess/remesh.hpp` — `repair_quality` (M-5)、界面適合 (M-4)。
- `src/oxidation.cpp` — Deal-Grove 定数。平面極限の較正に使う:
  線形速度定数 `k_s` と拡散係数 `D_ox` は B, B/A から
  `D_ox = B·N_ox/(2·C_gas)`、`k_s = (B/A)·N_ox/C_gas` で復元する
  (Deal-Grove 1965 の定義; N_ox = 2.25e22 cm⁻³ (SiO₂ 分子密度)、
  C_gas = 溶解酸化剤濃度 5.2e16 cm⁻³ (dry 1000°C) — 比だけが効くので
  C_gas は ParamDB `ox2d.cgas` 既定 5.2e16 とする)。

## 実装手順

### Stage A (コミット 1): 酸化剤定常拡散ソルバー

1. `include/cprocess/oxidant_solver.hpp` / `src/oxidant_solver.cpp` を新設:
   ```cpp
   struct OxidantResult {
     std::vector<double> conc;        // per-cell oxidant conc (oxide cells only)
     // per-face interface velocity: face index -> v_n [cm/s], Si 側へ正
     std::vector<std::pair<int,double>> iface_vn;
   };
   // Steady-state div(D grad C) = 0 in the oxide region.
   //   gas/oxide surface (oxide cell face on boundary patch zmax or facing a
   //   gas cell): Dirichlet C = c_gas
   //   oxide/Si interface face: Robin, flux = ks * C_face (consumed by reaction)
   //   oxide/nitride and all other faces: zero flux (nitride blocks oxidant)
   OxidantResult solve_oxidant(const Mesh& m, const std::vector<char>& oxide_mask,
                               const std::vector<char>& si_mask,
                               double d_ox, double ks, double c_gas,
                               std::ostream* log);
   ```
2. 実装は DiffusionSolver と同じ CSR + `cg_ilu0` パターンの縮小版
   (定常なので時間項なし):
   - 行 = 全セル。oxide_mask==0 のセルは恒等行 (C=0)。
   - 内部 oxide-oxide 面: 透過率 `tf = d_ox · g` (g は DiffusionSolver::build
     と同じ `|S|²/(S·d)` + 0.05 クランプ) で通常組立。
   - **Dirichlet (gas/oxide)**: `tb = d_ox · gb` を対角へ、`tb·c_gas` を rhs へ
     (DiffusionSolver::assemble の kBoundOwner 分岐と同形)。
   - **Robin (oxide/Si 界面面)**: 面フラックス `F = ks·C_face`、
     `C_face = C_P · gb·d_ox / (gb·d_ox + ks·|S|)`… ではなく標準の
     直列合成で組む: 実効透過率
     `t_robin = 1.0 / (1.0/(d_ox·gb) + 1.0/(ks·area))` (area = |S|) を
     対角に加算 (外側 C=0 への流出として扱う)。界面面の反応フラックス
     [個/s] は `F_f = t_robin · C_P`。**この式をコメントに明記。**
   - 解けたら各界面面で `v_n_f = F_f / (area · N_ox)` [cm/s] を返す。
   - 定常方程式は SPD なので `cg_ilu0` (rtol 1e-10, maxit 2000)。
3. **平面極限の検証テスト** (Stage A の DoD、tests/test_oxidant.cpp):
   マスクなし平板 (oxide 層厚 x_o 一様) で解いた v_n が Deal-Grove の
   `dx/dt = (B/2)/(x_o + A/2) / 60 / 60` (per-hour→per-s 換算) と
   相対差 < 10% で一致すること (x_o = 0.01 / 0.05 / 0.2 µm の 3 点)。

### Stage B (コミット 2): 幾何移動 (ノード毎 ALE) と oxidize への結線

4. `ale_move` のノード毎速度オーバーロードを追加
   (include/cprocess/ale_mover.hpp):
   ```cpp
   AleResult ale_move(Mesh& m, const MeshTopology& topo,
                      const std::vector<double>& vn_node,  // per-node, cm/s
                      double dt, std::function<bool(int)> mask_fn = nullptr,
                      int smooth_iters = 3);
   ```
   実装は現行関数の `v_n * dt` を `vn_node[i] * dt` に替えるだけ。
   法線は内部で `compute_node_normals` を呼ぶ (現行版もこの形に揃えて
   よいが既存シグネチャは残す)。
5. 面速度 → ノード速度の変換: 界面ノード i の `vn_node[i]` =
   そのノードを含む界面面の `v_n_f` の面積重み平均。非界面ノードは 0。
6. **`proc::oxidize` の 2D 分岐**: メッシュに `"nitride"` セルが存在する
   場合、P1-6/P2-3 の 1D 経路の代わりに次のループ (サブステップ N=10 の
   各ステップで):
   - (a) `solve_oxidant` (d_ox/ks は Deal-Grove から温度換算、上記 Stage A)
   - (b) **Si 消費**: 界面ノードを `vn_node` で Si 側 (−normal) に
     `ale_move(…, dt_k)` — 界面法線は Si→oxide 向きに正規化し、界面を
     Si 側へ `v_n·dt` 進める。
   - (c) **体積膨張**: 外面 (oxide/gas 境界) ノードを、対応する直下の
     界面ノード変位の **1.27 倍** だけ +z (外向き法線) に動かす。
     ノード対応は「同じ (x,y) 格子列」ではなく **ノード毎簿記** で行う:
     界面ノード j の変位 `d_j = v_n(j)·dt` を、外面ノード k のうち
     xy 距離最小の j に紐付け `rise_k = 1.27 · d_j` とする (最近傍探索、
     セル高より遠い場合は 0)。**この簿記をコメントで明記。**
     メッシュ上端に余裕が無い場合は先に `extend_mesh_exact` で
     ヘッドルーム 2 セル層を確保する。
   - (d) 界面通過セルの再タグ: 移動後、centroid が界面より oxide 側に
     入った silicon セルを ox_tag に (M-4 の界面適合を呼んでから)。
     nitride セルは決して再タグしない。
   - (e) `Mesh::finalize()` → 全フィールドポインタつき
     `repair_quality(st.mesh, &field_ptrs, 0.1)` → 体積変化の
     `rescale_fields_for_volume_change`。
   - (f) P2-3 の I 注入 (`theta·d_dx·N_Si/h`) は界面**面ごと**の
     `v_n_f·dt` を使ってセル単位に加算。点欠陥緩和 (diffuse_ted) も
     P2-3 と同様に dt_k で実行。
7. マスクなしの場合は従来 1D 経路のまま (回帰なし)。

## テスト仕様

`tests/test_oxidant.cpp` (Stage A) と `tests/test_ox2d.cpp` (Stage B) を
CMakeLists foreach に登録 (`oxidant`, `ox2d`)。

1. **平面 Deal-Grove 一致** (Stage A、上記手順 3): 3 膜厚点で v_n 相対差
   < 10%。さらに Stage B で: マスクなし + 2D 経路を強制 (テスト用に
   ダミー nitride を x>0.15µm に置き、開口部中央で測定) した
   `oxidize(30 min, 1000 °C, wet)` の開口中央酸化膜厚が
   `deal_grove_step` の解析値と相対差 < 15% (離散化込み)。
2. **bird's beak**: 0.4×0.2×0.5 µm メッシュ (8×4×50)、x∈[0.2,0.4] µm に
   nitride ストライプ (deposit("nitride", 0.05) + etch で x<0.2 を開口)、
   `oxidize(60 min, 1000 °C, wet)`。x 位置ごとの酸化膜厚 t(x)
   (oxide セルの z 範囲) を評価し:
   - マスク端から 0.1 µm 内側 (x = 0.3 µm) の t が開口部 (x = 0.1 µm) の
     **20〜80%**
   - t(x) がマスク下で単調非増加 (x = 0.22, 0.26, 0.30, 0.34 µm の 4 点で
     t[i+1] <= t[i]·1.05 — 5% のメッシュノイズ許容)
3. **メッシュ健全性**: テスト 2 の 5 サブステップ (N=10 のうちログで確認)
   を通して各ステップ後 `mesh_quality(st.mesh).min_q > 0.02`、反転セルなし
   (finalize が throw しない)。
4. **保存則**: テスト 2 で B 1e18 背景の総量変化 < 5%、全体積が
   「初期 + 1.27×消費 Si 体積 − 消費 Si 体積」と相対差 < 10%。

Python (`python/test_comprehensive.py`): deposit("nitride") + etch 開口 →
`oxidize(30, 1000, wet=True)` が完走し、開口部の表面 z が nitride 下より
高い (numpy で centroid フィルタ)。既存 oxidize テスト (1D 経路) PASS。

## 完了条件 (DoD)

- [ ] Stage A: `solve_oxidant` (Dirichlet/Robin/遮断 BC、上記式) + 平面極限
      テスト PASS — コミット 1 (`P2-4a`)
- [ ] Stage B: ノード毎 `ale_move` オーバーロード、1.27x 膨張簿記、
      oxidize の 2D 分岐、テスト 2〜4 PASS — コミット 2 (`P2-4b`)
- [ ] マスクなしフローは 1D 経路のままで既存テスト全 PASS
- [ ] Python テスト追加、`python3 python/test_comprehensive.py` PASS
- [ ] コミットメッセージに `P2-4` (a/b) を含める

## やらないこと

- **粘性流 / 応力による酸化遅延** (P3 (e); P2-6 の応力場が前提)
- Massoud 薄膜補正の 2D 化、面方位依存レート
- 窒化膜の酸化 (nitride は不透過・不変)
- narrow-band 化や levelset 表現 (界面は ALE ノード追跡; levelset は P2-5)
- メッシュ異方性細分化 (M-6)
