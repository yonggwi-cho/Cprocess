# M-3: 粗視化 (エッジ縮約)

## 目的

M-2 の適応細分化はセル数を単調増加させる。長時間アニール後に勾配が緩んだ領域の
セルを**エッジ縮約 (edge collapse)** で統合し、セル数の際限ない増大を抑える。
内部ライブラリ (remesh) のみの追加であり、proc:: 公開は将来の自動リメッシュ
ループ (P2) に委ねる。

## 前提知識 (実装者向け)

- **エッジ縮約 (b→a)**: エッジ (a,b) について、ノード b を a に併合する。
  a と b の**両方**を含むセル (エッジの incident セル) は体積ゼロに潰れるので
  削除。b のみを含むセルは接続の b を a に置換する。内部エッジのみを縮約し、
  境界/界面ノードを動かさない限り、置換セル群の体積和 = 削除前の体積和が
  厳密に成り立つ (メッシュ全体の総体積は保存される)
- **有効性検査 (簡易リンク条件)**: 完全なリンク条件の代わりに、
  「置換後の全セルの signed volume と品質を直接検査する」実用版で足りる
  (反転・退化・低品質をすべて弾けるため)

## 現状コード

- `src/remesh.cpp` / `include/cprocess/remesh.hpp` — `split_edges`,
  `redistribute_field`, `laplacian_smooth`, `flip_repair`, `repair_quality`。
  エッジ縮約は未実装。無名 namespace に `kEdges[6][2]` と
  `signed_vol(a,b,c,d)` あり (流用する)
- `include/cprocess/topology.hpp` — `MeshTopology::build(m)`:
  `edge_cells` (edge_key → incident cells), `node_cells`,
  `node_boundary` (境界面上のノード), `node_interface` (複数 region に接する
  ノード)。縮約の保護判定に使う
- `Mesh::finalize()` — cells/nodes から faces 等を再構築。縮約で未参照になる
  ノードが残っても finalize は動くが、ノード配列は本タスクでは**圧縮しない**
  (孤立ノードを許容する — スコープ外参照)

## 実装手順

1. `include/cprocess/remesh.hpp` に追加:
   ```cpp
   struct CoarsenResult {
     int n_collapsed = 0;   // 実行された縮約数
     int n_rejected = 0;    // 保護/品質/反転で拒否された数
     int n_cells_removed = 0;
   };
   // edges の各 (a,b) について b を a へ縮約する (a は動かない)。
   // fields が非 null なら列挙された全セル場を保存的に転写する。
   // 実行後に m.finalize() を内部で呼ぶ。
   CoarsenResult coarsen(Mesh& m,
                         const std::vector<std::pair<int,int>>& edges,
                         std::vector<std::vector<double>*>* fields = nullptr);

   // 縮約候補の自動選択: 全内部面のうち相対濃度差
   // |C_o - C_n| < rel_grad_thresh * max(C_o, C_n, 1e-3*global_max)
   // を満たす「低勾配面」の owner セルの最短エッジを候補とする。
   // 両端とも非境界・非界面ノードのエッジのみ。短い順にソートし、
   // 同一セルに触れる 2 本目以降を除外 (split_edges と同じ競合回避)。
   // 返す本数は全エッジ数 * max_fraction を上限とする。
   std::vector<std::pair<int,int>> select_coarsen_edges(
       const Mesh& m, const std::vector<double>& conc,
       double rel_grad_thresh, double max_fraction = 0.1);
   ```
2. `coarsen` の実装 (src/remesh.cpp)。処理順:
   - (a) `MeshTopology topo; topo.build(m);`、`vol0 = m.cell_vol` を退避。
     全フィールドの総質量 `M0_f = Σ C_f·V` を記録 (デバッグ assert 用)
   - (b) `std::vector<char> cell_touched(nc, 0)` を用意。各エッジ (a,b) を
     順に処理し、以下のいずれかで **拒否** (`++n_rejected`, continue):
     1. `topo.node_boundary[a] || topo.node_boundary[b] ||
        topo.node_interface[a] || topo.node_interface[b]`
        (境界/界面ノードは一切縮約しない — 外形と界面形状の保護)
     2. incident セル (`topo.edge_incident(a,b)`) が空、または
        incident/影響セルのどれかが `cell_touched`
     3. 置換シミュレーション: `topo.node_cells[b]` の各セルについて、
        (a,b) 両方を含むセルは「削除予定」、b のみのセルは b→a 置換後の
        signed volume を計算。**1 つでも volume <= 0 または
        `tet_quality` < 0.05 なら拒否** (これが簡易リンク条件を兼ねる:
        位相的に不正な縮約は必ず反転セルを生む)
   - (c) 受理したら実行: 削除予定セルに削除マーク、b のみのセルの接続を
     b→a に置換。関与した全セル (削除・置換とも) を `cell_touched` に
     マーク。ノード b はそのまま残す (孤立ノード、参照ゼロ)
   - (d) **場の転写** (fields 非 null 時)。規則を厳密に定める:
     1. 置換セル (b→a で形が変わったセル) は体積が変わる。各置換セルについて
        `C_new = C_old * V_old / V_new` とし、**セル単体の質量を保存**する
        (V_old は vol0、V_new は置換後座標から signed_vol で計算)
     2. 各削除セルの質量 `C_del * V_del` は、その削除セルと**ノードを 3 個
        以上共有する生存セル** (置換後接続で判定。b は a と同一視する) に
        **等分配**する: 対象セル k 本それぞれに質量 m/k を加算、すなわち
        `C_k += (C_del * V_del / k) / V_k_new`。対象がゼロ本の場合
        (起こらないはずだが) はエッジ incident 外の置換セル全体に等分配
     3. この 2 規則により全フィールドの総質量は厳密に保存される
        (丸めのみ)
   - (e) 全エッジ処理後、削除マークセルを erase して cells/cell_region/
     フィールドを圧縮 (flip_repair の compaction と同形)。`m.finalize()`。
     内部エッジのみ縮約なので `total_volume()` は厳密保存 —
     デバッグビルド用に `assert` 相当のチェックを入れてよい
3. `select_coarsen_edges` の実装:
   - `global_max = max(conc)`; 各内部面で上式の低勾配判定 → owner セルの
     **最短**エッジを候補化 (a<b 正規化、set で重複排除)
   - 両端の node_boundary/node_interface を事前フィルタ
   - エッジ長昇順にソート → 先頭から取りつつ、incident セル集合が既取得
     エッジと重なるものをスキップ → `総エッジ本数(edge_cells.size()) *
     max_fraction` 本で打ち切り
4. テストを `tests/test_remesh.cpp` に追加 (新規実行ファイル不要、
   CMakeLists.txt 変更なし)。**proc:: 追加なしのため Python 側の変更は不要**
   (docs/tasks/README.md 共通規約)

## テスト仕様

`tests/test_remesh.cpp` に追加:

1. **一様場の縮約**: 6×6×6 箱メッシュ (1×1×1 µm)、一様場 1e15 を 1 本用意。
   `select_coarsen_edges(m, conc, /*thresh*/1e9, /*max_fraction*/0.2)`
   (閾値を巨大にして全内部面を低勾配扱い) → `coarsen(m, edges, &fields)`。
   合格基準:
   - `n_collapsed >= 1` かつ `n_cells_removed > 0` (セル数が減る)
   - `total_volume()` の相対差 **< 1e-12**
   - 全フィールド質量 Σ C·V の相対差 **< 1e-9**
   - `mesh_quality(m).min_q > 0.05`
2. **勾配保護の選択**: テスト M-2 と同じ Gaussian B プロファイル
   (Rp=0.1µm, dRp=0.02µm) を 6×6×24 メッシュに設定。
   `select_coarsen_edges(m, B, 0.1, 0.2)` が返すエッジの**両端ノードいずれも**
   深さ帯 Rp±3dRp 内に z 座標を持たないこと (高勾配帯を一切選ばない)
3. **境界エッジの拒否**: 4×4×4 箱メッシュで、両端とも `node_boundary` の
   エッジだけを列挙して `coarsen` に渡す → `n_collapsed == 0`、
   `n_rejected == 候補数`、メッシュ完全不変 (セル数・総体積とも)
4. **非一様場の質量保存**: テスト 2 のメッシュ・場で遠方帯のエッジ 20 本を
   縮約 → B 総質量の相対差 < 1e-9、高勾配帯 (Rp±2dRp) のセル値の最大値が
   縮約前の最大値と相対差 < 1e-6 (プロファイルのピークが汚れない)
5. 既存の remesh / ale / integration テストがすべて PASS

## 完了条件 (DoD)

- [ ] `coarsen` / `select_coarsen_edges` 実装 (境界/界面保護、反転・品質
      ガード、2 段の保存的場転写、finalize 内包)
- [ ] 上記 4 テスト追加、`ctest --test-dir build` 全 PASS
- [ ] コミットメッセージに `M-3` を含める

## やらないこと

- proc:: / pybind / Simulation への公開 (内部ライブラリのみ。自動リメッシュ
  ループで公開するのは将来タスク)
- 孤立ノードの圧縮・再番号付け (finalize は孤立ノードを許容する。
  ノード配列のコンパクションは必要になった時点で別タスク)
- 完全なリンク条件 (Dey らの位相検査)。体積・品質の直接検査で代替する
- 境界・界面**上**のエッジ縮約 (半縮約含む)
- refine との自動交互ループ、並列化 (PA-1)
