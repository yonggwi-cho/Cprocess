# M-1: エッジ分割の完全自動化 + 一括分割 API

## 目的

現行 `edge_split(Mesh&, int a, int b)` (src/remesh.cpp) は cells /
cell_region を更新するが **faces・cell_vol・cell_cent 等を再構築しない**
(呼び出し側が `Mesh::finalize()` を呼ぶ暗黙の契約)。また 1 エッジずつの
O(nc) 走査で、多数エッジの分割が O(k·nc) になる。適応リメッシュ (M-2) と
品質修復 (M-5) の基盤として、(1) finalize までを内包した安全な API と
(2) 複数エッジの一括分割を実装する。

## 現状コード

- `src/remesh.cpp` — `edge_split(Mesh& m, int a, int b)`:
  全セルを走査して辺 (a,b) を含むセル (inc) を列挙 → 中点ノード追加 →
  各 inc セルを 2 分割 (b→mid 置換で in-place、a→mid 置換で push_back)。
  `cell_region` は複製するが faces は触らない。戻り値は新ノード番号 (なければ -1)
- `include/cprocess/topology.hpp` — `MeshTopology::build(m)` が
  `edge_cells` (エッジ→セル) マップを持つ (64bit エッジキー)。
  全セル走査の代わりにこれを使える
- `Mesh::finalize()` (src/mesh.cpp) — cells/nodes から faces, cell_vol,
  cell_cent, 境界パッチを再構築する
- フィールド (`SimState::fields`) は**セル番号に紐づく**。セル分割で
  セル数が増えるため、分割後のフィールド転写規則が必要

## 実装手順

1. `include/cprocess/remesh.hpp` に一括分割 API を追加:
   ```cpp
   struct SplitResult {
     int n_split = 0;          // 実際に分割されたエッジ数
     int n_skipped = 0;        // 競合でスキップされたエッジ数
     std::vector<int> cell_parent;  // 新メッシュの各セル -> 旧セル番号
   };
   // edges: (node_a, node_b) のリスト。分割後に m.finalize() 済みで返す。
   // 同一パス内で「同じセルに触る」2 本目以降のエッジはスキップし
   // n_skipped に計上する (呼び出し側が次パスで再試行する設計)。
   SplitResult split_edges(Mesh& m,
                           const std::vector<std::pair<int,int>>& edges);
   ```
2. `split_edges` の実装 (src/remesh.cpp):
   - (a) `MeshTopology topo; topo.build(m);` で `edge_cells` を作る
     (全セル走査を置換)
   - (b) `std::vector<char> cell_touched(nc, 0)` を用意。各エッジについて
     incident セルのどれかが touched ならスキップ (`++n_skipped`)。
     そうでなければ分割を実行し、incident セル (元セルと新セル両方) を
     touched にする
   - (c) 分割本体は現行 `edge_split` のループを流用。ただし
     `cell_parent` を維持する: 初期値 `iota(0..nc-1)`、in-place 更新セルは
     そのまま、push_back された新セルは `cell_parent.push_back(親の値)`
   - (d) 全エッジ処理後に `m.finalize()` を 1 回呼ぶ
   - (e) `SplitResult` を返す
3. 既存の単発 `edge_split(m, a, b)` は**後方互換のため残す**が、実装を
   `split_edges(m, {{a,b}})` への委譲に置き換える。戻り値仕様 (新ノード番号
   or -1) を維持するため、`SplitResult` に `last_new_node` を追加するか、
   委譲前後の `m.nodes.size()` 差分で判定する
4. フィールド転写ヘルパーを追加:
   ```cpp
   // conc (旧セル数) を cell_parent に従って新セル数へ複製転写する。
   // 濃度 [cm^-3] は intensive なので単純コピーで質量が保存される
   // (子セルの体積和 = 親セルの体積)。
   std::vector<double> redistribute_field(const std::vector<double>& conc,
                                          const std::vector<int>& cell_parent);
   ```
5. `tests/test_remesh.cpp` に追記 (新ファイル不要)

## テスト仕様

`tests/test_remesh.cpp` に追加:

1. **finalize 内包**: 4×4×4 箱メッシュで 1 エッジを `split_edges` →
   戻り後に `m.faces.size() > 0` かつ `m.cell_vol.size() == m.cells.size()`
   かつ `mesh_quality(m).min_q > 0` (反転なし)
2. **体積保存**: 分割前後で `m.total_volume()` の相対差 < 1e-12
3. **一括分割と競合スキップ**: 同一セルに属する 2 本のエッジを渡す →
   `n_split == 1 && n_skipped == 1`
4. **フィールド保存**: 一様濃度 1e15 を `redistribute_field` で転写 →
   全セル 1e15 のまま。さらに質量 Σ conc·vol が分割前後で相対差 < 1e-12
5. **多数エッジ**: 品質 > 0.1 の初期メッシュで 50 本をランダム選択
   (固定シード) して分割 → finalize 成功、min_q > 0、体積保存 < 1e-12
6. 既存の edge_split テストが (委譲実装でも) そのまま PASS

## 完了条件 (DoD)

- [ ] `split_edges` / `redistribute_field` 実装、`edge_split` は委譲に置換
- [ ] 上記 6 テスト追加、全テスト PASS
- [ ] `Mesh::finalize()` 呼び出しが `split_edges` 内に移動している
      (呼び出し側の手動 finalize が不要)
- [ ] コミットメッセージに `M-1` を含める

## やらないこと

- エッジ縮約・フリップ (M-3, M-4)
- 分割エッジの自動選択ロジック (M-2 適応細分化のスコープ)
- SimState/proc:: レベルの API (M-2 で `proc::refine` として公開する)
- 並列化 (PA-1 のスコープ。ここは正しさ優先で逐次)
