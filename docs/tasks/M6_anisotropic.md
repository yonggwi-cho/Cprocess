# M-6: 異方性メッシュ (方向依存細分化、長期)

## 目的

プロセスシミュレーションの解像度要求は本質的に異方的である: 濃度・界面は
深さ (z) 方向に急峻で、横方向には緩やか。等方細分化 (M-2) は縦解像度を
上げるたびに横方向のセルも巻き添えで増やし、セル数が爆発する。
本タスクは**方向依存の細分化** — 指定方向に整列したエッジのみを分割する —
を実装し、横方向の増加なしに層状の縦解像度を得る。
メトリックテンソルベースのフル異方性リメッシュは**スコープ外** (P3 水準の
別タスク) とし、実用的な directional refinement に限定する。

## 前提

- **M-2 (`refine_gradient` / `proc::refine`) と M-3 (`coarsen`) がマージ済み
  であること。** 本仕様は M-2 のパス構造・場転写・成長上限をそのまま流用する
- P3 ホライズンのタスクであり、着手時に M-2 実装の実際のシグネチャと本仕様の
  差分を確認してから開始する

## 現状コード (M-2/M-3 マージ後の想定)

- `src/remesh.cpp` — `refine_gradient(m, fields, key_index, thresh,
  max_passes, max_growth)`: 面走査 → 急勾配面の owner 最長エッジ分割 →
  `split_edges` → `redistribute_field` のパスループ
- `src/box_mesh.cpp` — `make_box_mesh`: Kuhn 6-tet 分割。箱メッシュでは
  エッジは軸方向 (x/y/z) と対角の混在
- `include/cprocess/process.hpp` — `proc::refine(st, species, thresh,
  passes, log)`
- `Vec3` (include/cprocess/vec3.hpp) — `dot`, `norm` 相当あり

## 実装手順

1. `include/cprocess/remesh.hpp` に追加:
   ```cpp
   // 方向依存細分化。M-2 の refine_gradient と同じパス構造だが、
   // エッジ指標が「勾配」ではなく「方向整列 × 長さ」:
   //   score(edge) = |dot(e_hat, d_hat)| * len(edge)
   // (e_hat: エッジ単位ベクトル、d_hat: direction 正規化)。
   // score > align_thresh * L_ref のエッジのみ分割する。
   //   L_ref = 現在のメッシュの「direction 方向エッジ」の平均長
   //           (|dot(e_hat,d_hat)| > 0.9 のエッジで平均。該当ゼロなら全平均)
   //   align_thresh 既定 0.5
   // fields は M-2 と同じ cell_parent 転写。max_growth も同じ意味。
   RefineResult refine_anisotropic(Mesh& m, const Vec3& direction,
                                   std::vector<std::vector<double>*>* fields,
                                   double align_thresh, int max_passes,
                                   double max_growth = 4.0);
   ```
2. 実装 (src/remesh.cpp)。各パス:
   - (a) `MeshTopology topo; topo.build(m);` の `edge_cells` を走査して
     全エッジを列挙 (a<b キーをデコード)
   - (b) `d_hat = direction / |direction|` (`|direction| <= 0` は
     `std::invalid_argument`)。各エッジで
     `score = |dot((n_b-n_a)/len, d_hat)| * len` を計算
   - (c) L_ref を上記定義で計算し、`score > align_thresh * L_ref` の
     エッジを候補化。**score 降順** (長く整列したものから) にソートして
     `split_edges` へ渡す (競合スキップは split_edges 任せ)
   - (d) fields を `redistribute_field` 転写、成長上限・パス上限は
     M-2 と同一。候補ゼロで終了
   - 注: 中点分割で生じる子エッジは半分の長さになり、数パスで
     `score <= align_thresh * L_ref` に落ちて自然に停止する (L_ref も
     パスごとに再計算されて縮む — これが「層を薄くする」駆動になる。
     したがって max_passes が実質の深さ制御であることを doc コメントに明記)
3. `proc::refine` に軸引数を追加する (**シグネチャ変更、後方互換は
   既定引数で確保**):
   ```cpp
   // axis: "" (既定) = M-2 の勾配駆動。 "x"|"y"|"z" = 方向依存細分化
   // (勾配指標は使わず refine_anisotropic に委譲。species は場転写の
   //  存在確認にのみ使用)。それ以外は std::runtime_error。
   void refine(SimState& st, const std::string& species,
               double rel_grad_thresh, int max_passes,
               const std::string& axis = "",
               std::ostream* log = nullptr);
   ```
   axis="z" のとき `direction = {0,0,1}`、align_thresh には
   rel_grad_thresh をそのまま渡す (既定 0.5 が両者で妥当)。
4. pybind (python/_cprocess.cpp): `proc_refine` に
   `py::arg("axis") = std::string("")` を追加。
5. Simulation (python/cprocess/simulation.py): `refine` に
   `axis: str = ""` キーワード引数を追加し `_c.proc_refine(...)` に渡す。
   docstring に「axis='z' で縦方向の層状細分化」と記載。単位変換なし。
6. テスト追加 (`tests/test_remesh.cpp` + `python/test_comprehensive.py`)。
   CMakeLists.txt 変更不要 (既存 remesh 実行ファイルに追記)。

## テスト仕様

`tests/test_remesh.cpp` に追加:

1. **縦解像度の向上**: 6×6×6 箱メッシュ (1×1×1 µm)、一様場 1e15。
   `refine_anisotropic(m, {0,0,1}, &fields, 0.5, 2)` 実行後:
   - 「縦エッジ」(|dot(e_hat, z_hat)| > 0.9) の**平均長**が実行前の
     **0.55 倍以下** (≈ 半減)
   - 「横エッジ」(|dot(e_hat, z_hat)| < 0.1) の**本数**の増加が
     実行前の **+10% 未満**
   - 場 1e15 が全セルで保たれ、質量相対差 < 1e-9、
     `total_volume()` 相対差 < 1e-12、`min_q > 0`
2. **方向の直交性**: 同メッシュで direction=(1,0,0) → 今度は x エッジ平均長
   が 0.55 倍以下、z エッジ本数増 < 10%
3. **不正入力**: `direction = {0,0,0}` で `std::invalid_argument` 送出
4. **proc 経由**: SimState + `proc::refine(st, "B", 0.5, 2, "z")`
   (B 一様場) がテスト 1 と同じ縦/横基準を満たす。`axis="w"` は
   `std::runtime_error`
5. 既存 remesh / ale / integration テスト全 PASS

`python/test_comprehensive.py` に追加:

- `sim.mesh(1,1,1,6,6,6)` → `sim.init("B",1e15)` → `n0 = sim.n_cells` →
  `sim.refine("B", passes=2, axis="z")` → `sim.n_cells > n0`、
  `sim.dose("B")` 相対差 < 1e-9、チェーン戻り値が sim

## 完了条件 (DoD)

- [ ] `refine_anisotropic` 実装 (整列×長さ指標、L_ref 自動基準、
      場転写・成長上限は M-2 流用)
- [ ] `proc::refine` の axis 引数 + pybind + Simulation.refine(axis=) +
      Python テスト
- [ ] 上記 C++ 4 テスト + Python 1 テスト追加、全テスト PASS
      (Python は `python3 python/test_comprehensive.py`)
- [ ] コミットメッセージに `M-6` を含める

## やらないこと

- メトリックテンソル場に基づくフル異方性リメッシュ (Hessian 指標、
  異方性 Delaunay 等) — 本タスクは方向指定 1 本の実用版
- 界面法線の自動推定による direction 決定 (呼び出し側が axis を指定する)
- 異方性**粗視化** (M-3 の等方縮約のまま)
- 勾配指標と方向指標の合成 (axis 指定時は勾配を見ない)
- 並列化
