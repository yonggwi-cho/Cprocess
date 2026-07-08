# P3-h: デバイスシミュレータ連携 (VTU 点データ + meta.json エクスポート)

## 目的

プロセス結果 (構造 + ドーピング) をデバイスシミュレータ (Sentaurus Device 等)
へ渡せる**中立フォーマット 2 本立て**で出力する (P3 概要の確定事項):

1. 既存 VTU 出力の拡張 — セル濃度場を**ノード平均した点データ**
   (`<PointData>`) として追加する。多くのデバイス側メッシュ変換ツールは
   節点値を要求するため
2. JSON サイドカー `<name>.meta.json` — 領域→材料表、境界パッチ名、単位系、
   種リスト、net doping の符号規約 (**ND − NA**、`DopType::neutral` は除外)

`proc::export_device(st, path_prefix, log)` が `<prefix>.vtu` +
`<prefix>.meta.json` の 2 ファイルを書く。TDR 直接出力はプロプライエタリの
ため対象外。**P1-11 (シリアライズ基盤) はマージ済みで前提。**

## 現状コード

- `src/vtk_writer.cpp` — `write_vtu(path, mesh, scalars, int_scalars)` (L9):
  ASCII XML UnstructuredGrid。座標は **µm に変換して出力** (`to_um = 1e4`,
  L17)、セルは VTK_TETRA(10) のみ、データは `<CellData>` の Float64/Int32
  だけで **`<PointData>` は未対応** (L48-66)。読込機能なし。
- `src/process.cpp` — `proc::save` (L2115): st.fields 全部 + 各ドーパントの
  `<Sym>_active` (`active_concentration` = C_ss クランプ) + `NetDoping`
  (donor は +、acceptor は −、**`DopType::neutral` は L2131 で除外済み**、
  P2-8) + int の `Region` を CellData で write_vtu に渡す。net doping の
  符号規約 ND−NA はここで既に実装されている — export_device は同じ組立を
  流用する。
- `include/cprocess/deck.hpp` — `SimState`: `mesh` / `fields`
  (map<string, vector<double>>, セル値) / `region_material`
  (map<int,string>) / `bcs` / `last_temp` / レジスト stack 一式。
- `include/cprocess/mesh.hpp` — `Mesh`: `nodes` (Vec3, cm) / `cells`
  (array<int,4>) / `cell_region` / `region_names` (map<int,string>) /
  `patch_names` (vector<string>) / `cell_vol` / `cell_cent` (finalize() で
  構築)。box メッシュのパッチ名は `{"xmin","xmax","ymin","ymax","zmin",
  "zmax"}` (src/box_mesh.cpp:48)。gmsh 入力では Physical Surface 名。
- `src/materials.cpp` — `kDopants` テーブル: B/In = acceptor、P/As/Sb =
  donor、C/F/Ge = **neutral** (P2-8)。`find_dopant(sym)` で st.fields の
  キーからドーパント種を判別できる (非ドーパント場: "I","V","C311",
  "damage", "sxx".. 等は nullptr)。
- `src/state_io.cpp` — CPRC1 バイナリ save/load (P1-11)。本タスクでは
  触らないが、`region_material` / `patch_names` の列挙パターンの参照実装。
- `CMakeLists.txt` — cprocess_core のソースリスト (L16-40) と
  テストの `foreach(t ...)` (L76)。
- `python/_cprocess.cpp` / `python/cprocess/simulation.py` — proc_*
  バインディングと `Simulation.save(path)` (パス passthrough、単位変換なし)
  の定型。

## 実装手順

1. **write_vtu の点データ対応** (vtk_writer.hpp / vtk_writer.cpp)。
   後方互換のデフォルト引数で拡張:
   ```cpp
   void write_vtu(
       const std::string& path, const Mesh& mesh,
       const std::vector<std::pair<std::string, const std::vector<double>*>>& scalars,
       const std::vector<std::pair<std::string, const std::vector<int>*>>& int_scalars,
       const std::vector<std::pair<std::string, const std::vector<double>*>>&
           point_scalars = {});
   ```
   `point_scalars` が空でなければ `</Cells>` の直後・`<CellData>` の前に
   `<PointData>` ブロックを書く (各配列 Float64, ascii, 値数 = n_nodes,
   既存と同じ `%.6e` 書式)。空なら出力は現行と**バイト単位で同一**
   (proc::save の回帰を保証)。

2. **セル→ノード体積加重平均** (新規 `include/cprocess/device_export.hpp` +
   `src/device_export.cpp`、CMakeLists の cprocess_core ソースリストで
   `src/state_io.cpp` の次に追加):
   ```cpp
   namespace cp {
   // Volume-weighted cell->node averaging:
   //   val(n) = sum_{c: n in mesh.cells[c], c not excluded} cell_vol[c]*f[c]
   //          / sum cell_vol[c]
   // exclude[c] != 0 marks cells skipped (gas); a node adjacent only to
   // excluded cells gets 0.
   std::vector<double> cell_to_node(const Mesh& mesh,
                                    const std::vector<double>& f,
                                    const std::vector<char>& exclude);
   // Writes <prefix>.vtu (cell data identical to proc::save + node-averaged
   // point data) and <prefix>.meta.json. Throws std::runtime_error on I/O.
   void write_device(const SimState& st, const std::string& prefix);
   }
   ```
   平均規則を厳密に定める: ノード n の値は、n を頂点に含む全セル c
   (`mesh.cells[c]` の 4 頂点走査で隣接リストを 1 回構築、O(4·nc)) のうち
   **gas セル (`region_material == "gas"`、material_id() 準拠の
   case-insensitive) を除外**した集合での `Σ V_c·f_c / Σ V_c`。
   分母 0 (gas にしか接しないノード) は 0.0。材料界面 (Si/oxide) を跨ぐ
   平均は許容 (デバイス側の慣行どおり両側平均。界面二重値は対象外 —
   やらないこと参照)。

3. **write_device の出力内容**:
   - **CellData**: proc::save (process.cpp L2115-2146) と同一の組立 —
     st.fields 全部、各ドーパントの `<Sym>_active`、`NetDoping`
     (`DopType::neutral` 除外、ND−NA)、int の `Region`。組立コードは
     proc::save から関数 `namespace proc { namespace detail { } }` に
     切り出さず、**process.cpp 側の proc::export_device で組み立てて
     cp::write_device に渡す**…ではなく単純化する: `write_device` が
     `find_dopant` / `active_concentration` を直接使って save と同じ規則で
     組み立てる (save 側は不変。二重実装だがテスト 4 で符号一致を担保)。
   - **PointData**: ドーパント濃度場 (find_dopant が非 null の fields キー) と
     その `<Sym>_active`、および `NetDoping` の 3 群を `cell_to_node` で
     ノード化し、**CellData と同名**で出力 (VTK は Point/Cell の同名配列を
     区別できる; ParaView でもそのまま両方見える)。"I"/"V"/"C311"/"damage"/
     応力 6 成分などの非ドーパント場は点データ化**しない** (セルデータの
     ままエクスポートはされる)。
   - **meta.json**: 手書き ofstream で下記スキーマを固定キー順で出力。
     文字列は `"` と `\` のみエスケープ (材料名・パッチ名は内部生成の
     ASCII)。数値は `%.17g`。

4. **meta.json スキーマ (全文例 — この例が本仕様の規範)**。2 領域
   (silicon 基板 + oxide 膜)・3 種 (B / As / C) の場合、
   `out.meta.json` は次の形になる:
   ```json
   {
     "format": "cprocess-device",
     "version": 1,
     "generator": "cprocess proc::export_device (P3-h)",
     "files": {
       "vtu": "out.vtu"
     },
     "units": {
       "length": "um",
       "concentration": "cm^-3",
       "temperature": "K",
       "note": "VTU point coordinates are in micrometres (write_vtu to_um=1e4); all concentration arrays are cm^-3"
     },
     "mesh": {
       "n_nodes": 1183,
       "n_cells": 4320,
       "cell_type": "tetrahedron"
     },
     "regions": [
       { "tag": 1, "name": "substrate", "material": "silicon" },
       { "tag": 2, "name": "film",      "material": "oxide"   }
     ],
     "boundaries": ["xmin", "xmax", "ymin", "ymax", "zmin", "zmax"],
     "species": [
       { "symbol": "B",  "name": "boron",   "type": "acceptor", "in_net_doping": true  },
       { "symbol": "As", "name": "arsenic", "type": "donor",    "in_net_doping": true  },
       { "symbol": "C",  "name": "carbon",  "type": "neutral",  "in_net_doping": false }
     ],
     "net_doping": {
       "field": "NetDoping",
       "convention": "ND-NA",
       "positive_means": "n-type",
       "uses_active_concentration": true,
       "excluded_types": ["neutral"]
     },
     "fields": {
       "cell": ["As", "B", "C", "I", "As_active", "B_active", "C_active", "NetDoping", "Region"],
       "point": ["As", "B", "C", "As_active", "B_active", "C_active", "NetDoping"]
     },
     "last_temp_k": 1173.15
   }
   ```
   生成規則: `regions` は `mesh.region_names` ∪ `region_material` のタグ
   和集合を tag 昇順で (名前が無ければ `"region<tag>"`、材料が無ければ
   `"silicon"` — silicon_mask の untagged=Si 規約と一致させる)。
   `boundaries` は `mesh.patch_names` をそのままの順で。`species` は
   st.fields のキーのうち `find_dopant` が非 null のものをキー昇順
   (std::map の走査順) で、`type` は DopType を "donor"/"acceptor"/
   "neutral" に写像、`in_net_doping = (type != neutral)`。`fields.cell` /
   `fields.point` は実際に VTU に書いた配列名を書いた順で列挙。

5. **`proc::export_device`** (process.hpp / process.cpp — 薄い委譲、
   P1-11 の save_state と同じ流儀):
   ```cpp
   // P3-h: device-simulator export. Writes <prefix>.vtu (cell data as
   // proc::save + node-averaged dopant/active/NetDoping point data) and
   // <prefix>.meta.json (region/material table, boundary patch names, unit
   // system, species list, ND-NA sign convention). Throws if a photoresist
   // stack is present (strip first) or on I/O error.
   void export_device(SimState& st, const std::string& path_prefix,
                      std::ostream* log = nullptr);
   ```
   `need_mesh(st)` → `st.has_stack` なら throw ("strip before export") →
   `cp::write_device(st, path_prefix)` → ログ
   `[export_device] wrote <prefix>.vtu + <prefix>.meta.json (<n_cells> cells, <n_species> species)`。

6. **pybind11 バインディング** (python/_cprocess.cpp、proc_* 規約):
   ```cpp
   m.def("proc_export_device",
       [](SimState& st, const std::string& prefix) {
         std::ostringstream log;
         proc::export_device(st, prefix, &log);
         return log.str();
       },
       py::arg("state"), py::arg("path_prefix"));
   ```

7. **`Simulation.export_device`** (python/cprocess/simulation.py、
   `save` の隣。パス passthrough、単位変換なし):
   ```python
   def export_device(self, path_prefix: str) -> "Simulation":
       """Device-simulator export (P3-h): writes <prefix>.vtu (with
       node-averaged point-data doping) and <prefix>.meta.json
       (regions/materials, boundaries, units, species, ND-NA convention)."""
       self._emit(_c.proc_export_device(self._st, path_prefix))
       return self
   ```

8. **テスト** (次節): `tests/test_device_export.cpp` 新規、CMakeLists の
   `foreach(t ...)` に `device_export` を追加。Python テストを
   `python/test_comprehensive.py` に追加 (`import json`,
   `import xml.etree.ElementTree as ET`)。

## テスト仕様

`tests/test_device_export.cpp` (一時ファイルはカレント相対
`"build_test_export.*"` で作り末尾に `std::remove`):

1. **save 回帰 (ParaView 互換の担保 1)**: 同一 state で `proc::save` した
   .vtu に文字列 `"<PointData>"` が**含まれない** (既存出力不変)。既存の
   save 利用テスト (`test_integration`, `test_flow`, `test_oxidize_flow`)
   が無変更で PASS することが本タスクの回帰条件。
2. **点データの整形式 (ParaView 互換の担保 2)**: mesh_box 1×1×1 µm
   (6×6×12) → `init("B",1e15)` → `implant_gauss("As", 1e15, 30 keV)` →
   `diffuse` (1000 °C, 10 min) → `export_device(st, "build_test_export")`。
   .vtu をテキストとして読み:
   - `<PointData>` ブロックが 1 個、`</PointData>` と対で存在し
     `<CellData>` より前にある
   - PointData 内の各 `<DataArray ... Name="As">` 等 (As, B, As_active,
     B_active, NetDoping の 5 本) の値行数が **n_nodes に一致**
   - `<DataArray` の開始タグ数と `</DataArray>` の終了タグ数が全文書で
     一致 (タグ収支による整形式チェック; 完全な XML 検証は Python 側
     ET.parse で行う)
3. **平均規則**: (a) 一様場 — 全セル B=3e17 のとき全ノード値が 3e17
   (相対差 < 1e-12; 体積加重平均の恒等性)。(b) gas 除外 — `etch` の
   polygon 版で一部を gas 化した state で、gas にしか接しないノードの
   点値が 0、Si/gas 界面ノードの点値が Si 側セルのみの加重平均に一致
   (手計算セルで検証、相対差 < 1e-12)。
4. **NetDoping の符号 (donor 注入 + p 基板)**: `init("B", 1e15)` (p 基板)
   + As Gaussian (dose 1e15, Rp=0.1 µm → ピーク ≈ 2e20 > 1e15)。
   `last_temp` を 1273.15 K にして export:
   - 表面近傍 (z_top から 0.1 µm 以内) のノード点値 NetDoping > 0 (n 型)
   - 深部 (z_top − 0.5 µm より深い) のノード点値 NetDoping < 0 (p 型、
     ≈ −1e15)
   - セル値 NetDoping が proc::save の出力と全セル相対差 < 1e-12
     (save との二重実装の一致担保)
   - `st.fields["C"] = 1e20` を加えて再 export → NetDoping (cell/point) が
     bit 不変 (**DopType::neutral の除外**、P2-8)
5. **meta.json 適合**: 上記 state (silicon + deposit した oxide の 2 領域、
   B/As/C の 3 種) の meta.json を読み、部分文字列で
   `"format": "cprocess-device"`, `"version": 1`, `"convention": "ND-NA"`,
   `"excluded_types": ["neutral"]`, `"material": "oxide"`,
   `"symbol": "C",` と同行の `"in_net_doping": false`, パッチ名 6 個、
   `"length": "um"` を CHECK (厳密パースは Python 側)。
6. **エラー系**: `photo` 後の export_device が throw、`strip` 後は成功。
   書込不能パス (`"/nonexistent_dir/x"`) で throw。mesh 無し state で throw。

`python/test_comprehensive.py` に追加 (`tempfile` + `json` + `ET`):

- フロー (mesh → init B → implant As → deposit oxide → diffuse) 後
  `sim.export_device(prefix)` がチェーン (`ret is sim`) で完走
- `meta = json.load(open(prefix + ".meta.json"))` が**パース成功** (スキーマ
  適合の operationalize) し、キー集合 ⊇ {"format","version","files","units",
  "mesh","regions","boundaries","species","net_doping","fields",
  "last_temp_k"}、`meta["net_doping"]["convention"] == "ND-NA"`、
  `meta["units"]["length"] == "um"`、species に As (donor, in_net_doping
  True) が含まれる、`meta["mesh"]["n_cells"] == sim.n_cells`
- `tree = ET.parse(prefix + ".vtu")` が**例外なく成功** (整形式 XML)、
  `tree.getroot().tag == "VTKFile"`、`.//PointData/DataArray` の Name 属性
  集合が `meta["fields"]["point"]` と一致
- 後始末 `os.unlink` × 2

**既存テストの回帰**: `test_state_io`, `test_integration`, `test_flow`,
`test_oxidize_flow`, `test_ox2d`, `test_mechanics` (write_vtu 利用箇所) が
無変更で PASS。`python3 python/test_comprehensive.py` PASS。

## 完了条件 (DoD)

- [ ] `write_vtu` に `point_scalars` デフォルト引数を追加、空なら出力
      バイト一致 (既存呼び出し無変更)
- [ ] `device_export.hpp/cpp` (`cell_to_node` 体積加重 + gas 除外、
      `write_device`) を cprocess_core に追加
- [ ] `<prefix>.meta.json` が本仕様掲載の JSON 例のスキーマ (キー・型・
      規約文字列) どおりに出力される
- [ ] `proc::export_device` + pybind `proc_export_device` +
      `Simulation.export_device` + Python テスト (json パース + キー検査 +
      ET.parse による VTU 再読込) — CLAUDE.md 必須ラッパ規約の 3 点セット
- [ ] `tests/test_device_export.cpp` 6 テスト PASS、CMakeLists foreach に
      `device_export` 登録、既存テスト全 PASS
- [ ] `build/_cprocess*.so` を `python/cprocess/` にコピーし
      `python3 python/test_comprehensive.py` PASS
- [ ] コミットメッセージに `P3-h` を含める

## やらないこと

- **TDR (Sentaurus 独自フォーマット) 直接出力** — プロプライエタリのため
  対象外 (概要の確定事項)。HDF5/CGNS/Exodus 等の他フォーマットも同様
- VTU のバイナリ/appended エンコーディング・圧縮 (ASCII のまま)、
  VTU の読込 (エクスポート専用)
- 境界パッチの**幾何** (面リスト) 出力 — meta.json は名前のみ。デバイス側
  はバウンディングボックスの面名 (xmin..zmax) で再構成できる前提
- 材料界面での点データ二重値 (界面ノードを材料ごとに分離した multi-valued
  export)。界面ノードは両側加重平均 1 値のみ
- 非ドーパント場 ("I"/"V"/"C311"/"damage"/応力 sxx..sxz/クラスタ *_cl) の
  点データ化 (セルデータとしては従来どおり出力される)
- proc::save 側の変更・リファクタ (CellData 組立の共通関数化はしない —
  一致はテスト 4 の < 1e-12 比較で担保)
- meta.json の外部 JSON Schema ファイル生成・バリデータ導入 (Python の
  json.load + キー assert で十分)
- キャリア統計 (n/p, 抵抗率) や ni 補正済み量の出力 — NetDoping (活性濃度
  ベース ND−NA) まで
- テキストデッキコマンド (`export_device file=...`) の追加 (Python 経由で
  十分。必要なら別タスク)
