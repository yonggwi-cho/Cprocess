# P1-11: 状態保存/再開 (バイナリ) + 1D 深さプロファイル抽出

## 目的

現状の出力は VTK .vtu の書き出しのみで読み戻せず、長いプロセスフローを
途中から再開できない。またプロファイル比較 (SIMS 相当) のたびに numpy で
手作業の抽出が要る。本タスクで (a) SimState の独自バイナリ save/load、
(b) 指定 (x,y) カラムの 1D 深さプロファイル抽出 `proc::profile1d`、
を実装する (SProcess の TDR 保存 + 1D extract の最小版)。

## 現状コード

- `src/vtk_writer.cpp` — `write_vtu(path, mesh, scalars, ints)`:
  ASCII VTK 書き出しのみ (読込なし)。本タスクでは触らない
- `include/cprocess/deck.hpp` — `SimState`: mesh / has_mesh / fields
  (map<string, vector<double>>) / region_material (map<int,string>) /
  bcs (vector<DirichletBC>: species, patch, conc) / last_temp /
  レジスト stack 一式 (has_stack, stack, ...)
- `include/cprocess/mesh.hpp` — Mesh: nodes (Vec3), cells (array<int,4>),
  cell_region, region_names, patch_names, faces ほか。faces/cell_vol 等は
  `finalize()` で再構築可能 → **保存しない**
- `src/process.cpp` — `proc::save` (vtu)。`need_mesh` ヘルパ
- `python/cprocess/simulation.py` — `Simulation.save` (vtu)

## 実装手順

1. 新規 `include/cprocess/state_io.hpp` + `src/state_io.cpp`
   (CMakeLists.txt の cprocess_core に追加) — 低レベル実装。
   proc:: 層は薄い委譲にする。
2. **バイナリフォーマット** (リトルエンディアン固定。x86/ARM64 前提で
   ホスト表現を直書きし、先頭に BOM 的マジックで判別):

   | # | 内容 | 型 |
   |---|------|----|
   | 1 | マジック `"CPRC1"` | char[5] |
   | 2 | フォーマットバージョン = 1 | int64 |
   | 3 | n_nodes | int64 |
   | 4 | ノード座標 x,y,z × n_nodes | double × 3·n_nodes |
   | 5 | n_cells | int64 |
   | 6 | セル接続 4 ノード id × n_cells | int64 × 4·n_cells |
   | 7 | cell_region × n_cells | int64 × n_cells |
   | 8 | n_region_names、続けて各 (tag int64, 名前長 int64, 名前 bytes) | — |
   | 9 | n_patch_names、続けて各 (名前長 int64, 名前 bytes) | — |
   | 10 | n_region_material、各 (tag int64, 材料名長 int64, bytes) | — |
   | 11 | n_fields、各 (種名長 int64, bytes, 値 double × n_cells) | — |
   | 12 | n_bcs、各 (種名長 int64, bytes, patch int64, conc double) | — |
   | 13 | last_temp | double |

   文字列はヌル終端なし。int はすべて int64 に拡幅して書く。
   フィールド長が n_cells と異なる場合は save 時に n_cells に
   resize(0.0 埋め) してから書く。
3. state_io API:
   ```cpp
   namespace cp {
   void write_state(const SimState& st, const std::string& path);  // throw on I/O error
   void read_state(SimState& st, const std::string& path);         // throw on I/O/format error
   }
   ```
   - `write_state`: `std::ofstream(path, std::ios::binary)` が開けなければ
     `std::runtime_error("cannot write '<path>'")`
   - `read_state`: マジック/バージョン不一致は
     `std::runtime_error("not a CPRC1 state file")` /
     `"unsupported state version"`。読込後:
     `st.mesh.finalize()`、`st.has_mesh = true`、レジスト stack 系メンバを
     すべて初期値にリセット (`st.has_stack=false` ほか)、layer_stack が
     存在するブランチでは clear (フォーマットには**含めない**)
4. **proc:: API** (process.hpp / process.cpp):
   ```cpp
   // レジスト stack は保存対象外: has_stack なら throw ("strip before save")。
   void save_state(SimState& st, const std::string& path,
                   std::ostream* log = nullptr);
   void load_state(SimState& st, const std::string& path,
                   std::ostream* log = nullptr);
   // (x_cm, y_cm) を xy-bbox に含むセルを集め、centroid.z 昇順の
   // (z_cm, conc) 列を返す。species が無ければ throw。該当セルゼロなら空。
   std::vector<std::pair<double,double>> profile1d(
       const SimState& st, const std::string& species,
       double x_cm, double y_cm);
   ```
   - save_state: `need_mesh(st)` → stack ガード → `write_state`。
     ログ `[save_state] wrote <path> (<n_cells> cells, <k> fields)`
   - load_state: `read_state`。ログ `[load_state] <path>: <n_cells> cells`
   - profile1d の**カラム所属判定**を厳密に定める: セル ci について
     4 頂点の x の min/max、y の min/max (= セルの xy-bbox) を計算し、
     `xmin <= x_cm <= xmax && ymin <= y_cm <= ymax` なら所属。
     全セルを線形走査 (O(nc) で十分)。gas 材料のセル
     (`region_material == "gas"`) は除外。所属セルを
     `cell_cent.z` 昇順に sort し `(cell_cent[ci].z, conc[ci])` を積む。
     同一 z 帯に複数セル (tet 分割で常に複数) があるのは許容 —
     呼び出し側でそのまま散布/平均できる
5. **pybind** (python/_cprocess.cpp):
   ```cpp
   m.def("proc_save_state", [](SimState& st, const std::string& p) {
       std::ostringstream log; proc::save_state(st, p, &log); return log.str(); },
     py::arg("state"), py::arg("path"));
   m.def("proc_load_state", [](SimState& st, const std::string& p) {
       std::ostringstream log; proc::load_state(st, p, &log); return log.str(); },
     py::arg("state"), py::arg("path"));
   m.def("proc_profile1d", [](const SimState& st, const std::string& sp,
                              double x, double y) {
       return proc::profile1d(st, sp, x, y); },  // list[(z_cm, conc)]
     py::arg("state"), py::arg("species"), py::arg("x"), py::arg("y"));
   ```
6. **Simulation** (python/cprocess/simulation.py):
   ```python
   def save_state(self, path: str) -> "Simulation":
       """Save the full simulation state to a binary CPRC1 file."""
       self._emit(_c.proc_save_state(self._st, path))
       return self

   def load_state(self, path: str) -> "Simulation":
       """Restore a state previously written by save_state()."""
       self._emit(_c.proc_load_state(self._st, path))
       return self

   def profile(self, species: str, x: float, y: float):
       """Depth profile at column (x, y) [um].
       Returns (z_um, conc) numpy arrays sorted by depth."""
       pairs = _c.proc_profile1d(self._st, species, x * UM, y * UM)
       z = np.array([p[0] for p in pairs]) / UM
       c = np.array([p[1] for p in pairs])
       return z, c
   ```
   save_state/load_state はパス passthrough (単位なし)、profile のみ
   µm↔cm 変換。profile はクエリメソッドなのでタプルを返す (チェーンしない)。
7. テスト: 新規 `tests/test_state_io.cpp` を作成し CMakeLists.txt の
   foreach に `state_io` を追加。一時ファイルは
   `"build_test_state.cprc"` 等カレント相対で作りテスト末尾に
   `std::remove`。

## テスト仕様

`tests/test_state_io.cpp`:

1. **roundtrip bit 一致**: mesh_box 1×1×1 µm (6×6×12) →
   `init("B",1e15)` → `implant_gauss("P", 1e13, energy=50)` →
   `deposit("oxide", 0.05µm)` → `add_bc("B", zmax パッチ, 1e17)` →
   `save_state` → 新しい SimState に `load_state`。合格基準:
   - n_cells / n_nodes が等しい、ノード座標・セル接続・cell_region が
     **bit 一致** (memcmp 相当)
   - fields のキー集合が一致し、全値 **bit 一致** (`==` で比較、許容差ゼロ)
   - region_material / region_names / patch_names / bcs / last_temp 一致
2. **再開後の solve**: load した状態で `diffuse` (1273.15 K, 60 s) が
   例外なく完走し、元の状態で同じ diffuse を実行した結果と
   全セル相対差 **< 1e-12** (finalize 再構築の決定性検証)
3. **profile1d**: 6×6×24 メッシュ (深さ 0.8µm)、B に縦 Gaussian
   (Rp=0.1µm, dRp=0.02µm) を設定 → `profile1d(st, "B", 中央 x, 中央 y)`:
   - 戻り列が空でなく z が**単調非減少**
   - conc 最大の要素の z が `z_top − Rp` と **1 セル高さ
     (0.8µm/24) 以内**で一致
4. **エラー系**:
   - `load_state(st, "no_such_file.cprc")` が throw
   - 先頭 5 バイトが `"XXXXX"` のダミーファイルで throw
   - `photo` 後 (has_stack) の `save_state` が throw、`strip` 後は成功
   - `profile1d(st, "Xx", ...)` (未知種) が throw、
     ドメイン外 (x = bbox 外) は空ベクトル (throw しない)

`python/test_comprehensive.py` に追加 (`import tempfile, os`):

- `tempfile.NamedTemporaryFile(suffix=".cprc", delete=False)` のパスに
  `sim.save_state(path)` → 新規 `Simulation` で `load_state(path)` →
  `n_cells` 一致、`np.array_equal(sim.field("B"), sim2.field("B"))`
  (bit 一致) → `sim2.diffuse(time=1, temp=1000)` が完走 →
  `z, c = sim2.profile("B", 0.5, 0.5)` で `len(z) > 0` かつ
  `np.all(np.diff(z) >= 0)` → `os.unlink(path)`

## 完了条件 (DoD)

- [ ] `state_io.hpp/cpp` (CPRC1 v1 フォーマット、上表の順序どおり) 実装、
      CMakeLists.txt の cprocess_core に追加
- [ ] `proc::save_state / load_state / profile1d` + pybind 3 本 +
      `Simulation.save_state / load_state / profile` + Python テスト
- [ ] load 後に `finalize()` が呼ばれ、stack 系メンバがリセットされる
- [ ] `tests/test_state_io.cpp` を foreach に登録、上記 4 テスト +
      Python テスト追加、全テスト PASS
- [ ] コミットメッセージに `P1-11` を含める

## やらないこと

- レジスト stack の保存 (save 前に strip を強制。フォーマットにも含めない —
  将来入れる場合は version 2)
- VTU の読込、TDR/HDF5 等の外部フォーマット、圧縮、エンディアン変換
  (リトルエンディアン環境のみサポートと文書化)
- profile1d の補間・等間隔リサンプリング・複数種同時抽出
  (numpy 側で加工できる生 (z, C) 列まで)
- デッキコマンド (`save_state file=...`) の追加 (Python 経由で十分。
  必要なら別タスク)
- 前方互換 (新しいコードで古い version のみ読める。逆は throw)
