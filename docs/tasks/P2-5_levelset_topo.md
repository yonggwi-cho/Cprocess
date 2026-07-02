# P2-5: レベルセット・トポグラフィ (レートベースのエッチ/デポ)

## 目的

現行のエッチ/デポは幾何指定 (深さ・厚さ + ポリゴン) のセル再タグであり、
等方性エッチのアンダーカットや共形デポのサイドウォール被覆を表現できない。
既存のレベルセットツールキット (`levelset.cpp`, 現在プロセスフロー未結線)
をエッチ/デポに結線し、材料別レート・等方/異方性を持つ
トポグラフィエンジンを実装する。

## 現状コード

- `src/levelset.cpp` — `levelset_init` (±1 シード + Dijkstra 近似符号付
  距離)、`levelset_advect` (一次風上 HJ、1 明示ステップ、|∇φ| は
  隣接ノード差分の max — 拡散的で不正確)、`levelset_reinit` (Sussman)、
  `levelset_update_regions` (φ の符号平均でセル再タグ)
- `src/process.cpp` — `etch`/`deposit` (幾何版、本タスク後も残す)、
  `extend_mesh_exact` (上方への headroom 追加に流用)
- `include/cprocess/topology.hpp` — `MeshTopology::node_adj`
- 前提: M-2 (適応細分化) マージ済み — 界面近傍の解像度確保に使う

## 実装手順

1. `levelset_advect` の勾配を改善する (`src/levelset.cpp`):
   - 現行: `|∇φ| ≈ max_j |φ_i − φ_j| / |x_i − x_j|` (隣接 j)
   - 改善 (ノード風上 Godunov 近似): 前進速度 F_i > 0 (界面が外向き) の
     ノードでは
     `g_i = max( max_j( (φ_i − φ_j)/d_ij, 0 ) )` を上流側差分として使い、
     F_i < 0 では符号を反転した式を使う (式を仕様に明記):
     ```
     F>0: |∇φ|_i = max_j max( (φ_i − φ_j)/d_ij, 0 )
     F<0: |∇φ|_i = max_j max( (φ_j − φ_i)/d_ij, 0 )
     ```
   - CFL: `dt_sub = 0.5 * h_min / max|F|` (h_min = 最短エッジ長) で
     サブステップ分割。`levelset_advect(m, topo, phi, F_node, dt)` の
     オーバーロードを追加 (F をノード場として受ける)
2. `proc::etch_rate` を追加 (`process.hpp`/`process.cpp`):
   ```cpp
   // rates: material name -> etch rate [cm/s]。載っていない材料は 0。
   // isotropic=false は垂直 (−z) 成分のみ (異方性/RIE 近似)。
   void etch_rate(SimState& st,
                  const std::map<std::string, double>& rates,
                  double time_s, bool isotropic,
                  std::ostream* log = nullptr);
   ```
   アルゴリズム:
   - (a) 露出表面から φ を初期化: gas セル (または境界上面) を正、
     材料セルを負として `levelset_init`
   - (b) ノード速度場 F_i: ノード近傍セルの材料のレートの最大値
     (multi-material 界面では速い方が食う)。異方性時は
     `F_i *= max(0, -n_z_i)` (n_z は M-5 系の `compute_node_normals` の
     z 成分; 上向き表面のみ削れる)
   - (c) `dt_sub` で time_s を刻み、advect → 5 サブステップごとに
     `levelset_reinit`
   - (d) 終了後 `levelset_update_regions` で φ>0 のセルを gas 再タグ、
     フィールドをゼロ化 (既存 etch と同じ規則)
   - (e) `repair_quality` は不要 (ノードは動かない — セル再タグのみ)。
     解像度不足は事前の `proc::refine` (M-2) 推奨、とログに助言を出す
3. `proc::deposit_conformal` を追加:
   ```cpp
   void deposit_conformal(SimState& st, const std::string& material,
                          double thickness_cm, std::ostream* log = nullptr);
   ```
   - (a) `extend_mesh_exact` で thickness 分の headroom (gas タグ) を確保
     (既存 deposit の機構を流用、追加セルは gas)
   - (b) 露出表面から φ 初期化、F = 一様 (thickness/time は無次元化し
     F=1, t=thickness で移流)
   - (c) φ が負に転じた gas セルを material に再タグ
   - 等方移流なのでステップ被覆 (サイドウォール) が自然に出る
4. pybind: `proc_etch_rate(state, rates_dict, time, isotropic)`,
   `proc_deposit_conformal(state, material, thickness)`。
   Python: `Simulation.etch_rate(rates, time, isotropic=True)`
   (rates は {"silicon": 0.1} µm/min → cm/s 変換)、
   `Simulation.deposit_conformal(material, thickness)` (µm)
5. デッキ: 追加しない (Python 専用。デッキは幾何版のみ維持)

## テスト仕様

`tests/test_topo.cpp` 新設 (CMakeLists foreach 登録):

1. **垂直エッチの一致**: 8×8×8 箱 Si、異方性 rate=0.1µm/min × 2min
   → 除去深さ 0.2µm が幾何 etch(depth=0.2) と同一セル集合 (± 1 セル層)
2. **等方性アンダーカット**: 上面の半分をレジスト (photo/mask で作成
   できないため region 再タグで nitride キャップを作る) で覆い、
   等方性エッチ深さ d 実行 → マスク端直下の横方向後退量が
   d の 50–130% (等方性の定性確認)
3. **共形デポ**: 幾何 etch でステップ (深さ 0.2µm、半面) を作り、
   deposit_conformal("nitride", 0.05µm) → ステップ側壁の中間高さに
   隣接する nitride セルが存在する (座標条件で assert)
4. **保存則**: etch_rate 後の全ドーパント質量が「残存セルの元の質量」に
   一致 (< 1e-9)
5. 既存テスト全 PASS (幾何 etch/deposit は無変更)

Python (`python/test_comprehensive.py`): etch_rate スモーク
(n_cells の gas 増加、field finite) + deposit_conformal スモーク。

## 完了条件 (DoD)

- [ ] 改善風上勾配 + CFL サブステップの advect オーバーロード
- [ ] `etch_rate` / `deposit_conformal` が proc::/pybind/Simulation 三層完備
- [ ] 上記 5 C++ テスト + Python テスト PASS、既存テスト回帰なし
- [ ] コミットメッセージに `P2-5` を含める

## やらないこと

- ノード移動 (ALE) による滑らかな表面追従 — セル再タグ解像度で表現
  (メッシュ細分化は M-2 の責務)
- プラズマ/RIE の物理モデル (イオン角度分布、ローディング効果)
- narrow-band 最適化・WENO 高次化 (性能タスクとして分離)
- 幾何版 etch/deposit の置換 (併存させる)
