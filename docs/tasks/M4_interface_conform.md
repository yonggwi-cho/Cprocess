# M-4: 局所再接続 (2-3 / 3-2 フリップ) と界面保護

## 目的

M-5 の「平滑化 + 分割」では直らない sliver がある (例: 針状 tet はどのエッジを
分割しても細長さが残る)。四面体メッシュの標準的な第 3 の修復手段である
**面フリップ (2-3) とエッジフリップ (3-2)** を実装し、`repair_quality` の
修復オペレータに加える。酸化 (P1-6) で界面が大きく動いてもメッシュ品質を
維持できるようにする。領域界面 (Si/SiO₂ 等) と境界パッチの面は**絶対に
フリップしない** (界面の形状が壊れるため) — これが「界面保護」の意味。

## 前提知識 (実装者向け)

- **2-3 フリップ**: 内部面 f を共有する 2 つの tet (頂点 {a,b,c,d} と
  {a,b,c,e}、共有面 {a,b,c}) を、辺 (d,e) を軸とする 3 つの tet
  {a,b,d,e}, {b,c,d,e}, {c,a,d,e} に置き換える。凸でない場合
  (線分 d-e が三角形 abc を貫かない場合) は新 tet に反転が出るので、
  **全新 tet の signed volume > 0 を検査して不可なら中止**すれば十分
- **3-2 フリップ**: 2-3 の逆。辺 (d,e) をちょうど 3 つの tet が共有する
  とき、それらを 2 つの tet に置き換える。実装は「辺 (d,e) の周囲セルが
  ちょうど 3 個」を `MeshTopology::edge_cells` で判定し、逆手順で構成
- どちらも適用条件は「置換後の全 tet の signed volume > ε」かつ
  「置換後の最小品質 > 置換前の最小品質」(品質改善フリップのみ実行)

## 現状コード

- `src/remesh.cpp` — `tet_quality`, `mesh_quality`, `laplacian_smooth`,
  M-1 の `split_edges`, M-5 の `repair_quality`
- `include/cprocess/topology.hpp` — `MeshTopology`: `edge_cells`
  (64bit キー edge→incident cells), `node_cells`, `node_boundary`,
  `node_interface`。face→cell は `Mesh::faces` (owner/neigh, neigh<0 が境界)
- `Mesh::cell_region` — 領域タグ。owner と neigh の region が異なる面が
  「界面」

## 実装手順

1. `include/cprocess/remesh.hpp` に追加:
   ```cpp
   struct FlipResult { int n_flip23 = 0, n_flip32 = 0; };
   // 品質 q < q_thresh のセルに接する面/辺に限定して 2-3 / 3-2 フリップを
   // 試みる。界面 (owner/neigh の region が異なる面) と境界面 (neigh < 0)
   // は候補から除外。フリップは「最小品質が改善する」場合のみ実行。
   // セル数が変わるため実行後に m.finalize() を内部で呼ぶ。
   // 注意: フリップはセルの合併/分割であり cell_parent での field 転写が
   // 定義できないため、out パラメータで「関与した旧セル集合と新セル集合」
   // を返し、呼び出し側 (repair_quality) が体積加重平均で転写する。
   struct FlipRemap { std::vector<int> old_cells, new_cells; };
   FlipResult flip_repair(Mesh& m, double q_thresh,
                          std::vector<FlipRemap>* remaps = nullptr);
   ```
2. 実装 (src/remesh.cpp)。処理順:
   - (a) `MeshTopology topo; topo.build(m);`
   - (b) sliver セル (q < q_thresh) を品質昇順に列挙
   - (c) 各 sliver セルについて:
     - **3-2 を先に試す** (セル数が減り品質改善が大きいことが多い):
       セルの 6 辺それぞれについて `edge_cells` の incident 数が 3 か確認 →
       3 なら 3-2 の置換 tet 2 個を構成し、(i) 全 signed vol > 1e-30、
       (ii) min(new q) > min(old q)、(iii) 3 セルすべて同一 region、
       (iv) 軸辺の 2 ノードが `node_boundary` でない、を満たせば実行
     - 3-2 不可なら **2-3**: セルの 4 面のうち内部面 (neigh >= 0) かつ
       owner/neigh 同一 region の面について、同様の 4 条件 (置換 3 tet の
       体積・品質改善・同一 region・面の 3 ノードすべて非境界) で実行
   - (d) フリップの実装はセル配列の直接編集で行う:
     置換前セルのうち 1〜2 個を新 tet で上書きし、余るセルは
     「削除マーク」(後で erase) または不足分を push_back。
     `cell_region` は関与セルの region (全て同一と保証済み) を引き継ぐ
   - (e) 1 パス終了後に削除マークのセルを消し、`m.finalize()`。
     フリップが 1 回でも起きたら topo は無効なので再 build するか
     1 パスで終える (**1 パスで終える実装を推奨**。繰り返しは
     repair_quality の round が担う)
3. `repair_quality` (M-5) の各 round に組み込む:
   平滑化 → 分割 → **flip_repair** → 平滑化 の順。FlipRemap を使った
   フィールド転写は repair_quality 内で行い、既存の `cell_parent` 合成と
   直交させる (フリップ関与セルは体積加重平均値を新セルに与える)
4. `RepairResult` に `n_flips` を追加

## テスト仕様

`tests/test_remesh.cpp` に追加:

1. **2-3 フリップ単体**: 2 つの tet {a,b,c,d}, {a,b,c,e} だけから成る
   手組みメッシュ (d, e が面 abc を挟んで両側、凸配置) を作り finalize →
   `flip_repair(m, 1.0)` (全セルが候補になる閾値) → `m.cells.size() == 3`、
   全セル体積 > 0、総体積が分割前と相対差 < 1e-12
2. **3-2 フリップ単体**: 上記の逆配置 (辺 de を 3 tet が共有) を手組み →
   `flip_repair` → `m.cells.size() == 2`、体積保存 < 1e-12
3. **界面保護**: テスト 1 の 2-tet 構成で `cell_region = {0, 1}` (異領域)
   に設定 → `flip_repair` が **n_flip23 == 0** (フリップしない)
4. **境界保護**: 箱メッシュ表面の面 (全 3 ノードが node_boundary) が
   フリップされないこと — 2×2×2 箱メッシュに q_thresh=1.0 で flip_repair を
   かけても総体積と外形 bbox が不変
5. **品質改善**: M-5 のテスト 1 と同じ人工 sliver メッシュで、
   `repair_quality` (フリップ組込み後) の `min_q_after` が
   フリップなし実装時の記録値以上
6. **フィールド保存**: フリップが起きるケースで一様濃度 1e15 →
   転写後も全セル 1e15、質量保存 < 1e-12

## 完了条件 (DoD)

- [ ] `flip_repair` 実装 (1 パス、界面/境界保護、品質改善フリップのみ)
- [ ] `repair_quality` に組込み、`RepairResult.n_flips` 報告
- [ ] 上記 6 テスト追加、全テスト PASS (既存 remesh/ale/integration 含む)
- [ ] コミットメッセージに `M-4` を含める

## やらないこと

- 4-4 フリップ・エッジ縮約 (M-3)
- Delaunay 性の保証 (品質改善のみを基準とする)
- 界面**上**の面のリメッシュ (界面は保護対象。界面形状の変更は
  ALE/levelset の仕事)
- 並列化 (PA-1)
