# PA-1: 残存逐次部の並列化 (リメッシュ候補収集 + MeshTopology::build)

## 目的

S-1 (ILU 分解並列化) 後に残る主要な逐次ホットスポットは、(a) リメッシュの
候補エッジ収集 (全セルの品質評価 + 最長エッジ走査)、(b) `MeshTopology::build`
の unordered_map 構築である。どちらもメッシュサイズに線形で、適応リメッシュ
(M-2) が繰り返し呼ぶため累積コストが大きい。両者を決定的 (スレッド数非依存)
に OpenMP 並列化する。内部ライブラリ変更のため Python 表面は不要。

## 現状コード

- `src/remesh.cpp` — `repair_quality(Mesh&, fields*, q_thresh, max_rounds,
  smooth_iters)`: ラウンドごとに全セルを逐次走査し、`tet_quality < q_thresh`
  のセルの最長エッジを求めて `std::set<std::pair<int,int>> seen` で重複排除
  しながら `candidates` に積む (~30 行)。この走査が O(nc) 逐次
- `src/topology.cpp` — `MeshTopology::build(const Mesh&)`: 逐次ループで
  (1) セル → `node_cells[c[v]].push_back(ci)` と
  `edge_cells[edge_key(a,b)].push_back(ci)` (unordered_map)、
  (2) edge_cells の全キーから `node_adj`、(3) 境界/界面フラグ。
  unordered_map への挿入が並列化を阻む
- `include/cprocess/topology.hpp` — `edge_cells` は
  `unordered_map<uint64_t, vector<int>>`。**検索はキー引き
  (`edge_incident`) のみで、反復順序に依存する利用箇所は node_adj 構築だけ**
- OpenMP パターンの手本: `src/mc_implant.cpp` (チャンク + スレッドローカル
  バッファ + 固定順逐次リダクション)、`src/sparse.cpp`

## 実装手順

### (a) repair_quality の候補エッジ収集の並列化

1. セル走査をスレッドローカル収集 + 決定的マージに置き換える:
   ```cpp
   const int nth = /* omp_get_max_threads(), 無効時 1 */;
   std::vector<std::vector<std::pair<int,int>>> tl(nth);
   #pragma omp parallel for schedule(static)
   for (int ci = 0; ci < nc; ++ci) {
     // 現行の品質評価 + 最長エッジ選択をそのまま (a<b 正規化まで)
     tl[omp_get_thread_num()].push_back({a, b});
   }
   ```
   `schedule(static)` により各スレッドの担当セル範囲は連続なので、
   `tl[0], tl[1], ...` の順に連結すればセル番号昇順が保たれる…とは
   **限らない** (chunk 境界はそうだが安全側で明示ソートする)。
   マージ後に `std::sort(candidates.begin(), candidates.end())` +
   `std::unique` で重複排除・順序確定する。これで **seen (std::set) は
   不要になり削除**、結果はスレッド数に依らず「セル番号順走査 + set」と
   同一集合・決定的順序 (辞書順) になる。
2. **注意**: 現行の `candidates` の並び順は「最初に発見したセル順」で、
   ソート後は「エッジ辞書順」に変わる。`split_edges` は先勝ちの競合スキップ
   をするため **分割結果が現行と変わり得る**。これは許容する (どちらも正当
   な決定的結果) が、シリアル版も同じソート順にすることで
   **1 スレッドと 4 スレッドの結果は厳密一致**させる。つまり順序規約は
   「エッジ (a,b) の辞書順」に統一する。
3. `tet_quality` 評価と最長エッジ選択はセルローカルなのでデータ競合なし。
   `omp.h` は `#ifdef _OPENMP` + スタブ (mc_implant.cpp の 3 行スタブ踏襲)。

### (b) MeshTopology::build の並列化

4. 2 パス方式に変更する:
   - **パス 1 (並列)**: スレッドローカルに
     `std::vector<std::pair<std::uint64_t,int>> tl_edges[nth]` を作り、
     `#pragma omp parallel for schedule(static)` のセルループで
     6 エッジ分の `(edge_key, ci)` を push する。`node_cells` は
     後段で作るのでここでは触らない
   - **パス 2 (逐次マージ)**: `tl_edges[0..nth-1]` をこの固定順で走査し
     `edge_cells[key].push_back(ci)` する。`schedule(static)` + スレッド順
     マージにより、各キーの cell リストは**セル番号昇順** (シリアルと同一)
     になる。map の**反復順序**は挿入順に依存し得るが、node_adj 構築後に
     `node_adj[i]` を `std::sort` することで順序依存を消す (現行は
     ソートしていない場合のみ追加; 呼び出し側はキー引きのみなので安全)
   - `node_cells`: ノード範囲分割の並列 2 パスで構築する。
     パス 1 (並列, セルループ) は各セルの 4 ノードをスレッドローカル
     `vector<pair<int,int>>` (node, cell) に積む。パス 2 はスレッド固定順の
     逐次 push (edge_cells と同じパターン) — cell 昇順が保たれる。
     もしくはより単純に「counts を並列で数えて offset を作り、逐次で
     詰める」でもよい (どちらでも決定的なら可)
   - `node_boundary` / `node_interface` ループはそれぞれ独立書込み
     (boundary はフラグ立てのみ → `#pragma omp parallel for` 可、
     interface はノードごと独立 → 並列 for) に変更
5. `node_adj` の構築は edge_cells 反復のままでよいが、末尾で各
   `node_adj[i]` をソートして順序を正規化する (手順 4 参照)。
6. 決定性の規約 (README 共通規約): 本タスクの並列化はすべて
   「スレッドローカル収集 → スレッド番号固定順マージ (or ソート)」であり、
   浮動小数を含まないため **1/4 スレッドで bit 一致** が要求水準。

## テスト仕様

`tests/test_remesh.cpp` と新規 `tests/test_topology_parallel.cpp` に追加。
CMakeLists の foreach に `topology_parallel` を追加。
`#ifdef _OPENMP` ガード + `omp_set_num_threads` はスタブパターン踏襲。

1. **topology 決定性**: `make_box_mesh(0,1, 0,1, 0,1, 24,24,24)` (24³) で
   `MeshTopology t1` (1 スレッド) / `t4` (4 スレッド) を build し、
   - エッジ総数 `t1.edge_cells.size() == t4.edge_cells.size()`
   - 全キーについて `t1.edge_cells[k] == t4.edge_cells[k]` (ベクトル完全一致)
   - 全ノードについて `node_cells` / `node_adj` / `node_boundary` /
     `node_interface` が完全一致
2. **repair_quality 決定性**: tests/test_remesh.cpp の M-5 sliver シナリオ
   (既存の repair_quality テストのメッシュ; 無ければ box メッシュの 1 ノードを
   面近くへ動かして sliver を作る) を 1 スレッドと 4 スレッドで実行し、
   `RepairResult` の `n_split` / `n_flips` / `n_smoothed` が一致、かつ
   最終メッシュのセル数・ノード数が一致すること。
3. **リグレッション**: 既存の test_remesh / test_mesh / test_ale /
   test_integration が PASS (候補順序の規約変更で分割位置が変わった場合、
   品質基準ベースの assert は満たされるはず。座標を直接 assert している
   テストがあれば期待値を更新し、その旨をコミットメッセージに書く)。
4. **時間 (情報のみ)**: 24³ の build を 1/4 スレッドで `std::chrono` 計測し
   printf。assert なし。

## 完了条件 (DoD)

- [ ] 候補エッジ収集と MeshTopology::build が OpenMP 並列化され、
      OpenMP 無効ビルドでも動く
- [ ] 1 スレッド vs 4 スレッドで bit 一致 (テスト 1・2)
- [ ] 既存全テスト + 新テスト PASS、CMakeLists 登録
- [ ] コミットメッセージに `PA-1` を含める

## やらないこと

- `split_edges` 本体・フリップ本体の並列化 (競合判定が逐次依存;
  彩色ベース並列フリップは別タスク)
- `edge_cells` のデータ構造変更 (ソート済みフラット配列化等) —
  効果はあるが API 互換を壊すため別途
- laplacian_smooth の並列化 (逐次スイープの意味論が変わるため)
- Python バインディング
