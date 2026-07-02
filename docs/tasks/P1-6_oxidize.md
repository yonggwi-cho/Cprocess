# P1-6: 酸化のプロセスフロー統合 (`proc::oxidize`)

## 目的

酸化は現在 `src/oxidation.cpp` の Deal-Grove **厚さ計算式のみ**で、メッシュも
SiO₂ 領域も動かない (`examples/nmos_locos_flow.cpp` が手動で結線したデモが
あるだけ)。本タスクで `proc::oxidize(time, temp, wet)` を新設し、
「露出した Si 上面の 1D ブランケット酸化」をデッキ / Python から 1 コマンドで
実行できるようにする。内容: 既存酸化膜厚の計測 → Deal-Grove 増分 →
体積収支 (SiO₂ は消費 Si の約 2.2 倍; 表面は 0.56·Δx_ox 上昇、Si は
0.44·Δx_ox 消費) → `extend_mesh_exact` + セル再タグで形状実現 →
全フィールド転写 → `repair_quality` → `st.last_temp` 更新。

## 現状コード

- `src/oxidation.cpp` / `include/cprocess/oxidation.hpp` —
  `deal_grove_step(x0_um, dt_min, T_celsius, wet)`。**単位に注意**:
  引数は µm / 分 / **摂氏**、戻り値は新しい酸化膜厚 (µm、増分ではなく総厚)。
  内部で dt を時間 (hr) に換算する (レート定数が per-hour のため)。
  コアの cm / s / K から呼ぶときは
  `x_um = x_cm * 1e4`, `dt_min = time_s / 60.0`, `T_c = temp_k - 273.15`
  と換算する。
- `src/process.cpp` — 無名 namespace の `extend_mesh_exact(base, thickness,
  nz_add)`: 元の nx/ny/nz を `infer_box_dims` で復元し、上方に `thickness` を
  足した箱メッシュを `make_box_mesh` で作り直す (z 方向は均等再分割)。
  `deposit()` がこのパターンの完全な手本: 合成領域タグ
  (`max_tag + 1000`) の付与、centroid 判定による再タグ、最近傍 centroid に
  よる全フィールド転写、`st.region_material` 更新。
- `src/process.cpp` — `photo()` / `infer_box_dims` / `point_in_polygon` /
  `silicon_mask` / `is_silicon`。`etch()` は gas タグに `max_tag + 2000` を
  使っている (衝突回避のため oxidize は `+3000` を使う)。
- `include/cprocess/remesh.hpp` — `repair_quality(Mesh&,
  std::vector<std::vector<double>*>* fields, q_thresh=0.1, max_rounds=3,
  smooth_iters=5)` (フィールド転写つきオーバーロード、M-5/M-4 で実装済み)。
- `include/cprocess/field_transfer.hpp` — `transfer_field_nearest`,
  `check_mass_conservation`。
- `include/cprocess/ale_mover.hpp` — 本タスクでは**使わない** (下記
  「やらないこと」参照)。
- `examples/nmos_locos_flow.cpp` — 手動結線デモ (Step 5–7)。今回の製品化対象。
  本タスク完了後もデモはそのまま残す (変更不要)。
- `src/deck.cpp` — `cmd_diffuse` が `Unit::time` / `Unit::temp` 換算の手本。
  `run_deck` のコマンドディスパッチに追加する。
- `include/cprocess/deck.hpp` — `SimState::last_temp` (K)。
- Python: `python/_cprocess.cpp` の `proc_deposit` バインディング、
  `python/cprocess/simulation.py` の `deposit()` メソッドが手本
  (`UM = 1e-4`, `MIN = 60.0`, `_celsius_to_k`)。

## 実装手順

1. `include/cprocess/process.hpp` に宣言を追加:
   ```cpp
   // Blanket 1D vertical thermal oxidation of the exposed top surface.
   // time_s in seconds, temp_k in K. Grows/extends the SiO2 layer above the
   // silicon per Deal-Grove; consumes 0.44*dx_ox of Si and raises the outer
   // surface by 0.56*dx_ox. Returns the new total oxide thickness in cm.
   double oxidize(SimState& st, double time_s, double temp_k, bool wet = false,
                  std::ostream* log = nullptr);
   ```
2. `src/process.cpp` に実装 (アルゴリズム、すべて cm 単位):
   - (a) `need_mesh(st)`; `time_s <= 0` は throw。
   - (b) **既存酸化膜厚の計測**: `z_si_top` = 材料が silicon
     (`is_silicon`、未タグ含む) のセル全体での「セルの 4 ノード z の最大値」の
     最大。`z_top = st.mesh.bbox().hi.z`。`z_si_top` より上に centroid を持つ
     セルを走査し、材料が `"oxide"`/`"sio2"` 以外 (gas 等) のセルが 1 つでも
     あれば `throw std::runtime_error("oxidize: top surface is not bare Si or
     SiO2")`。既存膜厚 `x0_cm = z_top - z_si_top` (酸化膜セルが無ければ 0;
     数値誤差対策に `x0_cm < 1e-9` なら 0 とする)。
   - (c) **Deal-Grove 増分**:
     ```cpp
     const double x_new_um = deal_grove_step(x0_cm * 1e4, time_s / 60.0,
                                             temp_k - 273.15, wet);
     const double dx_ox = x_new_um * 1e-4 - x0_cm;   // grown oxide, cm
     ```
     `dx_ox <= 0` なら形状変更なしでログのみ出して return。
   - (d) **体積収支と形状実現**: 消費 Si 厚 `dx_si = 0.44 * dx_ox`、表面上昇
     `rise = 0.56 * dx_ox`。z セル高
     `h = base_lz / nz` (`infer_box_dims` の nz、`base_lz = bbox` の z 幅)。
     `nz_add = std::max(1, (int)std::round(rise / h))` で
     `Mesh ext = extend_mesh_exact(st.mesh, rise, nz_add);` を呼ぶ。
   - (e) **再タグ**: 新 Si/SiO₂ 界面高さ `z_if = z_si_top - dx_si`。
     酸化膜タグは、`st.region_material` に既に `"oxide"` の合成タグ
     (>= 1000) があればそれを再利用、無ければ `max_tag + 3000` を新設して
     `st.region_material[ox_tag] = "oxide"`。`ext.cell_region` を deposit()
     と同様に centroid で再構成する:
     - `c.z > z_if` → `ox_tag` (旧表面より上の新規酸化膜 + Si から転化した
       上部 0.44·dx_ox 帯の両方を含む)
     - それ以外 → 旧メッシュの最近傍 centroid セルの `cell_region` をコピー
       (deposit() の転写ループと同じ最近傍探索を region にも使う)。
   - (f) **フィールド転写**: deposit() と同じ最近傍 centroid 転写で
     `st.fields` の全フィールドを `ext` に移す。その後、**旧表面より上**
     (`c.z > z_si_top + x0_cm` すなわち旧 bbox 上面より上) の新規セルは
     ドーパントを 0 にする。**Si→SiO₂ に転化したセル
     (`z_if < c.z <= z_si_top`) は転写値をそのまま保持する** — この帯の
     ドーパントは P1-4 の偏析で扱う (それまでは非 Si マスクで凍結される)。
   - (g) `st.mesh = std::move(ext);`
   - (h) **品質修復**: `st.fields` の全フィールドへのポインタ配列を作り
     `RepairResult rr = repair_quality(st.mesh, &field_ptrs, 0.1);` を呼ぶ
     (extend_mesh_exact は正則箱メッシュを返すので通常 no-op で即 return
     するが、将来の非正則メッシュに備え必ず呼ぶ)。`rr.min_q_after` をログ。
   - (i) `st.last_temp = temp_k;`
   - (j) ログ例:
     `[oxidize] dry 1273.15 K 1800 s: tox 0 -> 0.054 um (dSi=0.0238 um,
     rise=0.0302 um), mesh now N tets, min_q=...`。膜厚は必ず
     `fmt("%.4g", x_new_um)` 形式で出す (Python テストが参照)。
   - 戻り値: `x_new_um * 1e-4` (cm)。
3. **デッキコマンド** (`src/deck.cpp`):
   ```
   oxidize time=30min temp=1000C ambient=dry|wet
   ```
   `cmd_oxidize`: `time` は `Unit::time`、`temp` は `Unit::temp`、
   `ambient` は省略可 (既定 dry)。`"dry"`/`"wet"` 以外は `c.fail(...)`。
   `run_deck` のディスパッチ連鎖に `else if (c.name == "oxidize")` を追加。
4. **pybind11 バインディング** (`python/_cprocess.cpp`、CLAUDE.md 規約):
   ```cpp
   m.def("proc_oxidize",
       [](SimState& st, double time_s, double temp_k, bool wet) {
         std::ostringstream log;
         proc::oxidize(st, time_s, temp_k, wet, &log);
         return log.str();
       },
       py::arg("state"), py::arg("time_s"), py::arg("temp_k"),
       py::arg("wet") = false,
       "Blanket thermal oxidation (Deal-Grove). time in s, temp in K.\n"
       "Returns a log string.");
   ```
5. **Simulation メソッド** (`python/cprocess/simulation.py`、deposit() の直後):
   ```python
   def oxidize(self, time: float, temp: float, *, wet: bool = False) -> "Simulation":
       """Blanket thermal oxidation of the exposed silicon top surface.

       time in minutes, temp in Celsius. Grows SiO2 per Deal-Grove
       (<100> Si); the surface rises by 0.56x and silicon is consumed
       by 0.44x of the grown oxide thickness.
       """
       self._emit(_c.proc_oxidize(self._st, time * MIN,
                                  _celsius_to_k(temp), bool(wet)))
       return self
   ```
6. `CMakeLists.txt` の `foreach(t ...)` に `oxidize_flow` を追加。

## テスト仕様

`tests/test_oxidize_flow.cpp` (新規)。共通メッシュ: 0.2×0.2×0.5 µm、
`mesh_box(st, 0, 0.2e-4, 0, 0.2e-4, 0, 0.5e-4, 4, 4, 50)` (z セル高 0.01 µm)。

1. **Deal-Grove 解析値との一致 (dry 1000 °C 60 min)**:
   `proc::oxidize(st, 3600.0, 1273.15, false)` の戻り値 tox が
   解析値 **0.0540 µm** (= `deal_grove_step(0, 60, 1000, false)`;
   B(1273 K)≈1.043e-2 µm²/hr, A≈0.139 µm) と相対差 < 1e-9 で一致
   (同じ関数を呼ぶだけなので厳密)。さらに**タグから測った**酸化膜厚
   (oxide セル centroid z の最大セルの上面 − 最小セルの下面; 実装は
   `bbox().hi.z − z_si_top` の再計測でよい) が 0.0540 µm の **±20% 以内**
   (セル高 0.01 µm の離散化で ±1 層 ≈ ±18% を許容)。
2. **繰り返しの自己整合性**: 新しい state で
   `oxidize(1800 s)` を 2 回 ↔ 別 state で `oxidize(3600 s)` を 1 回。
   両者の戻り値 tox の相対差 < 5% (Deal-Grove の厳密積分なので実際は
   ~1e-12 だが、既存膜厚計測がメッシュ離散化を経由するため 5% とする)。
3. **体積収支**: oxidize 前後の全セル体積和の増分が
   `0.56 * dx_ox * (0.2e-4)^2` (= rise × 底面積) と相対差 < 5%
   (nz_add の丸めが唯一の誤差源)。
4. **ドーパント保存 (Si→SiO₂ 帯の凍結保持)**: oxidize 前に
   `proc::init(st, "B", 1e18)`。oxidize 後、(a) 全メッシュの B 総量
   (Σ C·V) が前後で相対差 < 5% (最近傍転写の smearing 許容)、(b) 旧表面
   (0.5 µm) より上のセルの B は全て 0、(c) `z_if < z <= 0.5 µm` の
   oxide 転化セルには B が残っている (max > 0.5e18)。
5. **異常系**: `etch` で表面を gas にした state で oxidize →
   `std::runtime_error` が投げられる。

デッキ: `tests/test_oxidize_flow.cpp` 内で `run_deck` に
`"oxidize time=60min temp=1000C ambient=dry"` を含む小デッキを食わせ、
例外なく完走しログに `[oxidize]` が含まれること。

Python (`python/test_comprehensive.py` に追加):

1. `Simulation().mesh(0.2, 0.2, 0.5, 4, 4, 50).init("B", 1e18)
   .oxidize(60, 1000)` が通り、`sim.cell_centroids[:,2].max()` が
   0.5 µm 超 (表面上昇) かつ増分が 0.56×0.054 µm = 0.0302 µm の ±50%
   以内 (centroid 判定のため緩め)。ログに `[oxidize]` を含む。
2. `oxidize(60, 1000, wet=True)` の膜厚 (ログの数値 or 2 回目の
   centroid 増分) が dry より大きい。
3. 既存の Python テストが全て PASS。

既存 C++ テスト (`ctest --test-dir build`) が全て PASS すること。

## 完了条件 (DoD)

- [ ] `proc::oxidize` 実装 (process.hpp/process.cpp、cm/s/K API、戻り値 = 総膜厚 cm)
- [ ] `deal_grove_step` の µm/min/°C 換算が上記どおり明示的に書かれている
- [ ] デッキコマンド `oxidize` (time/temp/ambient) 追加
- [ ] `proc_oxidize` バインディング + `Simulation.oxidize` (min/°C 変換、`self` 返し)
- [ ] `tests/test_oxidize_flow.cpp` の 5 テスト + デッキテストが PASS、
      `CMakeLists.txt` の foreach に `oxidize_flow` 追加
- [ ] Python テスト追加、`python3 python/test_comprehensive.py` PASS
- [ ] `cmake --build build -j$(nproc) && ctest --test-dir build` 全 PASS
- [ ] コミットメッセージに `P1-6` を含める

## やらないこと

- **ALE によるノード移動での滑らかな界面移動** (P2-4 の 2D/3D 酸化で実施)。
  本タスクは `extend_mesh_exact` + centroid 再タグの階段状実現で十分
- Massoud 薄膜補正、面方位/圧力/HCl 依存レート (P2-4)
- マスク開口つき局所酸化 (LOCOS / bird's beak) — ブランケットのみ (P2-4)
- Si/SiO₂ 界面の偏析・dose loss (P1-4; 本タスクでは転化セルの濃度を凍結保持
  するだけ)
- OED/ORD (P2-3)、酸化応力 (P2-6)
- `examples/nmos_locos_flow.cpp` の書き換え
