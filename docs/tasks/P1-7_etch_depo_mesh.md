# P1-7: エッチ/デポの表面適合改善 (多層スタック + 真の材料除去 + 選択比)

## 目的

現行の deposit/etch は「箱メッシュ再生成 + セル再タグ」であり、
(a) 多層デポでの層管理が region_material の合成タグ任せで壊れやすい、
(b) エッチは gas 再タグのままセルが残り、以後の bbox・注入表面 (z_top) が
物理表面とずれる、(c) 材料選択比がない。本タスクで
(a) `SimState::layer_stack` による層帳簿、(b) ブランケットエッチの
**真のセル除去** (メッシュ縮小)、(c) `etch(material=...)` 選択エッチ、
を実装する。表面の平滑化・レート物理は P2-5、酸化との層整合は P1-6
(依存として名前参照のみ) のスコープ。

## 現状コード

- `src/process.cpp`:
  - `infer_box_dims(m)` — ノード座標の一意値カウントで (nx,ny,nz) を復元
  - `extend_mesh_exact(base, thickness, nz_add)` — bbox から
    `make_box_mesh(..., z0, z_top_old + thickness, nx, ny, nz + nz_add)`
    を再生成。**bbox ベースなので 2 回目の deposit でも「現 top の上」に
    正しく積まれる** (実装確認事項 — テスト 1 で検証する)。ただし層の
    素性は `region_material[dep_tag]` の合成タグ (max_tag+1000) のみで、
    2 層目は max_tag+1000 が再計算されて別タグになる — 帳簿がない
  - `deposit(st, material, thickness, nz_add, poly, log)` — 上記 + 全場を
    最近傍 centroid で O(nc²) 転写 (この転写は本タスクで触らない)
  - `etch(st, depth, poly, log)` — centroid z > z_cut かつポリゴン内の
    セルを gas タグに再タグし濃度をゼロ化。**セルは残る**
- `include/cprocess/deck.hpp` — `SimState` (layer_stack なし)
- `python/_cprocess.cpp` — `proc_deposit` / `proc_etch`
- `python/cprocess/simulation.py` — `Simulation.deposit` / `Simulation.etch`
- 依存 (名前参照のみ): `docs/tasks/P1-6_oxidize.md` (酸化フロー統合)。
  P1-6 がまだ無い場合、oxidize 連携箇所は「maintained by deposit/etch
  (+ 将来 oxidize)」というコメントに留める

## 実装手順

1. **層帳簿** — `include/cprocess/deck.hpp` の `SimState` に追加:
   ```cpp
   // 上から順の (region_tag, material)。deposit が先頭に push、
   // ブランケットエッチが層を消費したら pop。mesh_box/mesh_gmsh で clear。
   // (将来 P1-6 oxidize も維持する)
   std::vector<std::pair<int,std::string>> layer_stack;
   ```
   - `mesh_box` / `mesh_gmsh` で `st.layer_stack.clear()`
   - `deposit` 成功時に `st.layer_stack.insert(begin(), {dep_tag, mat})`
   - `deposit` の既存バグ修正: `for (int tag : ext.region_tags())
     st.region_material[tag] = "silicon";` は**既存の非 Si 層のタグを
     silicon で上書きする**。修正: 既に `st.region_material` に存在する
     タグは上書きしない (存在しないタグのみ "silicon" を設定)。
     さらに ext の cell_region 再構成で、旧メッシュの層 (z <= z_top) は
     **最近傍 centroid の旧セルの region を引き継ぐ**必要がある —
     現状は 0 に初期化されたままなので、既存の場転写ループと同じ最近傍
     探索で `ext.cell_region[ci] = st.mesh.cell_region[best_j]` (z <= z_top
     のセルについて) を設定する
2. **etch のシグネチャ拡張** (process.hpp / process.cpp):
   ```cpp
   // material: 空文字なら全材料 (gas を除く) を対象。非空なら
   // その材料 (region_material の値が lower(material) に一致、
   // "si"/"silicon" 等のエイリアスは is_silicon と同様に正規化) の
   // セルのみ除去/再タグ対象とする。
   void etch(SimState& st, double depth,
             const std::vector<std::pair<double,double>>& poly = {},
             const std::string& material = "",
             std::ostream* log = nullptr);
   ```
3. **ブランケットエッチ = 真の除去** (poly が空 or <3 頂点のとき):
   - (a) `z_cut = bbox().hi.z - depth`。除去対象 =
     「centroid.z > z_cut」かつ「material フィルタに合致」のセル。
     部分深さセルの規則: **centroid が z_cut より上なら除去、下なら保持**
     (階段状は許容 — 平滑表面は P2-5)
   - (b) 新メッシュ構築: 保持セルのみで `Mesh` を作り直す。
     使用ノードを収集して**ノード圧縮** (old→new のインデックスマップ)、
     `cells` / `cell_region` を保持セル分だけコピー、`patch_names` は
     旧メッシュからコピー、`finalize()`。
     フィールド転写は**保持セルの恒等コピー** (`new_conc[k] =
     old_conc[kept[k]]`) — 補間なし、保持セルの質量は厳密保存
   - (c) 除去で消滅した層タグ (保持セルに 1 つも残らない region) を
     `st.layer_stack` から erase (region_material のエントリは残してよい)
   - (d) `st.has_stack` (レジスト) が true なら
     `std::runtime_error("etch: strip resist first")`
4. **ポリゴンエッチ = 従来どおり gas 再タグ** (poly 指定時):
   挙動は現行のまま (material フィルタのみ追加適用)。
   **除去にしない理由 (コメントに明記すること)**: 内部カラムだけをセル除去
   すると上面が非多様体 (穴あき/オーバーハング) になり、box-mesh 前提の
   `extend_mesh_exact` / `infer_box_dims` (次の deposit) が破綻するため。
   トポグラフィの一般化は P2-5 (レベルセット) のスコープ
5. **選択エッチ**: material 非空のとき、深さ帯 (z > z_cut, poly 内) の
   セルのうち **region_material が一致するものだけ**を除去 (ブランケット)
   または gas 再タグ (ポリゴン)。一致しないセルは濃度も含め不変。
   ブランケット選択エッチで「対象材料の上に別材料が残る」ケース
   (埋め込み層) は対象外 — 深さ帯判定のみで良い
6. **pybind** (python/_cprocess.cpp): `proc_etch` に
   `py::arg("material") = std::string("")` を追加 (poly の後)
7. **Simulation** (python/cprocess/simulation.py): `etch` に
   `material: str = ""` キーワードを追加し passthrough。docstring 更新
   (「poly なしのブランケットエッチはセルを物理的に除去し bbox が縮む。
   poly 指定時は gas 再タグ」)。deposit の docstring に多層可を追記
8. **デッキ** (src/deck.cpp): `etch` コマンドがあれば `material=` を追加、
   なければ変更不要 (deck.hpp のコマンド一覧に etch が無いことを確認済み
   なら何もしない)

## テスト仕様

新規 `tests/test_etch_depo.cpp` を作成し、CMakeLists.txt の
`foreach(t ...)` に `etch_depo` を追加:

1. **多層デポ**: mesh_box 1×1×0.5 µm (nx=ny=4, nz=8) →
   `deposit("oxide", 0.1µm)` → `deposit("nitride", 0.1µm)`。合格基準:
   - bbox.hi.z が初期 +0.2 µm (相対差 < 1e-9)
   - `layer_stack.size() == 2`、`layer_stack[0].second == "nitride"`、
     `layer_stack[1].second == "oxide"`
   - z 帯 (top−0.1, top) の全セルの region_material が "nitride"、
     z 帯 (top−0.2, top−0.1) が "oxide" (centroid で判定)
   - 1 層目 (oxide) のタグが 2 回目の deposit 後も region_material 上で
     "oxide" のまま (silicon 上書きバグの回帰テスト)
2. **ブランケットエッチの除去**: 上記の後 `etch(0.1µm)` →
   - `n_cells` が減少し、bbox.hi.z が (top − 0.1µm) と
     **1 セル高さ以内**で一致
   - `layer_stack.size() == 1` (nitride が消えて oxide が最上)
   - `mesh_quality(m).min_q > 0`、`cell_vol.size() == cells.size()`
3. **ドーパント質量保存**: mesh_box → `init("B", 1e15)` →
   `implant_gauss` (B, dose=1e13, energy=50) → 深いエッチ前に
   「z <= z_cut のセルの B 質量」を記録 → `etch(0.05µm)` →
   残メッシュの B 総質量が記録値と相対差 **< 1e-9** (恒等コピーなので
   実際は < 1e-12)
4. **選択エッチ**: mesh_box → `deposit("oxide", 0.1µm)` →
   Si セル数を記録 → `etch(0.2µm, {}, "oxide")` →
   - oxide セルが 0 になる (深さ帯 0.2µm > 膜厚 0.1µm)
   - **Si セル数が不変** (深さ帯にかかる Si セルも material 不一致で保持)
   - bbox.hi.z が元の Si 上面と 1 セル高さ以内で一致
5. **ポリゴンエッチは従来挙動**: `etch(0.1µm, 正方形ポリゴン)` →
   n_cells **不変**、ポリゴン内・深さ帯のセルが gas タグ、濃度ゼロ。
   既存 integration テスト (etch 使用箇所) がそのまま PASS
6. **レジストガード**: photo 後の blanket etch が throw

`python/test_comprehensive.py` に追加:

- `sim.deposit("oxide", 0.1).deposit("nitride", 0.1)` 後
  `sim.etch(0.1)` で `sim.n_cells` 減少・bbox 縮小、続けて
  `sim.etch(0.2, material="oxide")` で例外なく動作し `sim.n_cells` 減少。
  `sim.etch(0.05, poly=[(0.2,0.2),(0.8,0.2),(0.8,0.8),(0.2,0.8)])` は
  n_cells 不変 (再タグ)。dose("B") の保存 (< 1e-9) も 1 ケース

## 完了条件 (DoD)

- [ ] `SimState::layer_stack` 追加、deposit/etch が維持、mesh_* が clear
- [ ] deposit の region 上書き/引継ぎバグ修正 (実装手順 1)
- [ ] blanket etch のセル除去 (ノード圧縮 + 恒等場転写)、
      polygon etch は再タグ維持、material 選択比
- [ ] pybind `proc_etch(material=)` + `Simulation.etch(material=)` +
      Python テスト、docstring 更新
- [ ] `tests/test_etch_depo.cpp` を CMakeLists.txt の foreach に登録、
      上記 6 テスト + Python テスト追加、全テスト PASS
- [ ] コミットメッセージに `P1-7` を含める

## やらないこと

- ポリゴンエッチの真の除去 (非多様体表面になるため — 実装手順 4 の理由)
- 表面平滑化・テーパ・等方/異方レート・時間発展 (P2-5 レベルセット)
- deposit の共形成長 (現状の上方積層のみ)、O(nc²) 最近傍転写の高速化
- 酸化との層整合 (P1-6 側で layer_stack を維持する)
- 部分深さセルの体積分割 (centroid 判定の階段で確定)
