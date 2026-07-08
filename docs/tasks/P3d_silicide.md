# P3-d: シリサイド化 (`proc::silicide`)

## 目的

Ni/Ti 堆積 + 熱処理での NiSi/TiSi₂ 形成を 1 コマンド化する。成長則は
Deal-Grove 類似の拡散律速 x² = x₀² + B·Δt (B は金属種/相ごとの Arrhenius、
NiSi は 1.5 eV オーダー)。Si 消費 (NiSi: 金属 1 に対し Si 1.83 消費) の
体積 bookkeeping は oxidize と同じ再タグ方式、ドーパントのシリサイドへの
偏析 (Si 側の dose loss) は P1-4/P1-9 の偏析機構 (SegTable) に材料対と
係数を追加して実現する (P3_overview.md の確定方針)。

**P1-6 (oxidize の消費/成長機構) と P1-7 (layer_stack) の完了が前提
(両方マージ済み)。**

## 現状コード

- `src/materials.cpp` / `include/cprocess/materials.hpp` — **一般の
  「材料表」は存在しない**。材料は (1) `MatId` enum (`kMatSi=0 ..
  kMatGas=4`) + `material_id(name)` の名前マッピング (未知名は kMatGas)、
  (2) `st.region_material` の材料名文字列、(3) **Dopant 構造体の
  材料別カラム** (`dox0/eox`, `dnit0/enit`, `dpoly0/epoly` — P1-9) の
  3 点で表現される。よって「nickel/nisi を材料表に追加する」の実体は
  **MatId 追加 + material_id() 追加 + Dopant にシリサイド用カラム追加**
  である (下記手順 1)。
- `include/cprocess/diffusion.hpp` — `SegTable` は **`double h[5][5]` /
  `m[5][5]` の固定サイズ** (MatId で添字)。MatId を増やすには配列次元の
  拡張が必須。`src/diffusion.cpp` の `make_seg_table()` (L298) が
  `for (i<5)` で h を埋め、`t.m[kMatSi][kMatOxide] = segregation_m(dp,T)`
  だけを非 1 にしている — シリサイド対はここに 1 行足す。
- `src/diffusion.cpp` — `DiffusionSolver` コンストラクタ (L58):
  `cell_mat` が与えられると `mask_[i] = (mat_[i] != kMatGas)` で solve
  マスクを再導出する。**新 MatId (metal/silicide) は自動的に「活性」に
  なる**が、`assemble()` の kSegregation 分岐 (L439) に
  「D<=0 の側があればフラックス遮断」ガードがあるため、金属の D=0 は
  そのまま完全バリア、シリサイドは D>0 を与えたときだけ偏析交換が働く。
- `src/process.cpp` — `oxidize()` の P1-6 レガシーパス (L1210–1308) が
  「厚さ計測 → 増分 → 界面高さ z_if での centroid 再タグ」の手本。
  ただし酸化は表面が上昇するのに対し **シリサイド化は表面が沈む**
  (NiSi: 1−r_Si−r_met = 1−0.82−0.45 = −0.27 < 0)。よって本タスクは
  extend_mesh_exact を使わず **既存 bbox 内の再タグのみ** で実現できる
  (最上部の余り帯を "gas" に再タグ = 表面後退; polygon etch() の
  gas 再タグ方式 L616 以降と同じ前例)。メッシュ再構築なし =
  フィールド転写も不要 (gas 帯のゼロ化のみ)。
- `src/process.cpp` — 合成タグの慣例: photo/deposit `+1000`、etch
  `+2000`、oxidize `+3000`。**silicide は `max_tag + 4000` を使う**。
  `deposit()` の `known[]` 材料リスト (L510) と `set_region()` の
  `known[]` (L214) — 金属名を追加する。
- `src/process.cpp` — blanket etch の layer_stack 掃除
  (`tag_has_cells` で消えた層を erase、L700–709) — 金属層が完全消費
  されたときの掃除に同じパターンを使う。
- `src/oxidation.cpp` — `deal_grove_step`: 線形項 (A) 込みなので流用
  **しない**。シリサイドは純拡散律速 x=√(x₀²+B·t) を silicide() 内で
  直接計算する (2 行なので専用ソースファイルは作らない)。
- ParamDB (P1-10): `ParamDB::instance().get(key, fallback)` パターン。
  materials.hpp 冒頭のキー一覧コメントに新キーを追記する。
- Python: `python/_cprocess.cpp` の `proc_oxidize` バインディング、
  `simulation.py` の `oxidize()` (L292) が手本。

## 実装手順

1. **材料の追加** (`include/cprocess/materials.hpp` / `src/materials.cpp`):
   - `MatId` に追記 (既存値は不変、gas の後ろに追加):
     ```cpp
     enum MatId : int { kMatSi = 0, kMatOxide = 1, kMatNitride = 2,
                        kMatPoly = 3, kMatGas = 4,
                        kMatSilicide = 5,   // NiSi / TiSi2 (P3-d)
                        kMatMetal = 6 };    // nickel / titanium (P3-d)
     constexpr int kMatCount = 7;
     ```
     単一の kMatSilicide を両相で共有する (oxide が sio2 を包含するのと
     同じ; 相の区別は region_material の文字列 "nisi"/"tisi2" が保持)。
   - `material_id()`: `"nickel"/"ni"` と `"titanium"/"ti"` → kMatMetal、
     `"nisi"/"tisi2"` → kMatSilicide。未知名 → kMatGas は不変。
   - `SegTable` (diffusion.hpp) を `h[kMatCount][kMatCount]` /
     `m[kMatCount][kMatCount]` に拡張。m の全 1 初期化は 5x5 の
     手書きブレースをやめ、デフォルトコンストラクタのループで埋める。
     `make_seg_table()` のループ上限も kMatCount に。追加行:
     ```cpp
     t.m[kMatSi][kMatSilicide] = segregation_m_silicide(dp, temp_k);
     ```
     (Si/metal 対は m=1 のままだが金属 D=0 の per-face ガードで遮断)。
   - `Dopant` 構造体の**末尾**にシリサイド用カラムを追加 (末尾追加なら
     kDopants の位置指定集約初期化は無修正でデフォルト値が入る):
     ```cpp
     // ── Silicide (P3-d): NiSi/TiSi2 共通の実効値 ──
     // 偏析平衡比 m_sil(T) = C_si / C_silicide = seg_sil_m0*exp(-seg_sil_e/kT)
     // (< 1 = シリサイド側優先 -> Si の dose loss)
     double seg_sil_m0 = 0.3, seg_sil_e = 0.0;
     // シリサイド中拡散係数 D_sil = dsil0*exp(-esil/kT) [cm^2/s]
     // (0 だと per-face D>0 ガードで偏析交換ごと遮断されるので既定は非 0)
     double dsil0 = 1.0e-3, esil = 2.0;
     ```
   - アクセサ追加 (segregation_m と同型、ParamDB 経由):
     `segregation_m_silicide(d, T)` (キー `<Sym>.seg_sil_m0` /
     `<Sym>.seg_sil_e`)、`silicide_diffusivity(d, T)` (キー
     `<Sym>.dsil0` / `<Sym>.esil`)。materials.hpp 冒頭のキー一覧
     コメントに 4 キーを追記。
   - `material_diffusivity()` の switch に
     `case kMatSilicide: return silicide_diffusivity(d, temp_k);` と
     `case kMatMetal: return 0.0;` (完全バリア) を追加。
   - 回帰確認: `grep -n "\[5\]\[5\]\|< 5" src include` で MatId=5 固定の
     残置がないこと (SegTable/make_seg_table 以外にはない想定)。
2. **材料名の受理**: `deposit()` の `known[]` に `"nickel","ni",
   "titanium","ti"` を追加 (金属はユーザが deposit する)。`set_region()`
   の `known[]` にも同 4 名を追加。"nisi"/"tisi2" は silicide() だけが
   生成する (deposit 不可のまま)。
3. `include/cprocess/process.hpp` に宣言を追加 (oxidize_2d の直後):
   ```cpp
   // Blanket silicidation (P3-d): converts a previously deposited blanket
   // metal film ("nickel" or "titanium") on the exposed Si top surface into
   // its silicide ("nisi" / "tisi2") by diffusion-limited growth
   //   x^2 = x0^2 + B(T)*time_s,  B = b0*exp(-eb/kT)   (ParamDB, per phase).
   // Volume bookkeeping (all retag-only, no mesh rebuild -- the surface
   // *recedes* by (rsi+rmet-1)*dx): growing dx of silicide consumes
   // rsi*dx of Si (interface moves down) and rmet*dx of metal; the excess
   // band at the top is retagged "gas". Growth stops when the metal is
   // exhausted (x capped at x0 + t_metal/rmet). Returns the new total
   // silicide thickness in cm. metal: "nickel"|"ni"|"titanium"|"ti".
   double silicide(SimState& st, const std::string& metal, double temp_k,
                   double time_s, std::ostream* log = nullptr);
   ```
   (引数順 metal, temp, time は P3_overview の確定シグネチャ)。
4. `src/process.cpp` に実装 (cm/s/K、単発ステップ — oxidize の
   theta=0 レガシーパス相当。OED 型サブステップはしない):
   - (a) `need_mesh`; `st.has_stack` / `time_s <= 0` は throw。metal 名を
     lower して相を決定: nickel→("nisi", キー接頭辞 "silicide.nisi")、
     titanium→("tisi2", "silicide.tisi2")。他は throw。
   - (b) **層計測** (oxidize L1181 と同型): `z_si_top` = is_silicon セルの
     最大ノード z。z_si_top より上の centroid のセルを走査し、許容材料は
     {当該金属, 当該相のシリサイド, gas} のみ; 金属セルが 1 つも無い、
     または他材料 (oxide 等) が混在する場合は
     `throw std::runtime_error("silicide: top surface is not <metal> on Si")`。
     金属帯厚 `t_m` = 金属セルの (最大ノード z − 最小ノード z)、既存
     シリサイド厚 `x0` = 同相シリサイドセル帯厚 (無ければ 0、1e-9 未満は 0)。
   - (c) **成長則**: ParamDB から
     `b0 = get("<pfx>.b0", 既定)`, `eb = get("<pfx>.eb", 既定)`,
     `rsi = get("<pfx>.rsi", 既定)`, `rmet = get("<pfx>.rmet", 既定)`
     (既定値は下表)。`B = b0*exp(-eb/(kBoltzmannEv*temp_k))`、
     `x_new = sqrt(x0*x0 + B*time_s)`。**金属律速キャップ**:
     `x_new = std::min(x_new, x0 + t_m / rmet)` (キャップ発動をログに
     `metal fully consumed` と明記)。`dx = x_new - x0`; `dx <= 0` なら
     ログのみで return。
   - (d) **再タグ (メッシュ再構築なし)**: 増分の帯割り:
     消費 Si 厚 `d_si = rsi*dx`、消費金属厚 `d_met = rmet*dx`、表面後退
     `d_gas = (rsi + rmet - 1.0)*dx` (NiSi 0.27dx / TiSi₂ 0.30dx、両相で
     正)。新界面 `z_if = z_si_top - d_si`。シリサイドタグは
     region_material に既存の同相合成タグ (>=1000) があれば再利用、
     無ければ `max_tag + 4000` を新設し
     `st.region_material[sil_tag] = "nisi"|"tisi2"`。centroid z で再タグ:
     - `z_if < z <= z_if + x_new` → sil_tag
     - `z_if + x_new < z <= z_top_old - d_gas` → 金属タグ (残存金属)
     - `z > z_top_old - d_gas` → gas (etch(poly) と同じ gas 再タグ;
       既存の gas 合成タグ再利用 or `max_tag + 2000` 慣例に合わせ新設)
     - `z <= z_if` → 不変 (触らない)
   - (e) **フィールド処理**: メッシュ不変なので転写なし。gas に再タグ
     したセルの全フィールドを 0 化。**Si→silicide 転化セルのドーパントは
     保持** (凍結ではなく、手順 1 の SegTable/D_sil により後続 diffuse で
     偏析交換され dose loss が現れる — oxidize の P1-4 連携と同じ設計)。
   - (f) **layer_stack**: sil_tag が新設なら
     `st.layer_stack.insert(begin, {sil_tag, "nisi"|"tisi2"})` を金属層
     エントリの位置 (金属より下) に挿入; 金属が完全消費されタグにセルが
     残らない場合は blanket etch と同じ `tag_has_cells` 掃除で erase。
   - (g) `st.last_temp = temp_k;` 戻り値 `x_new` (cm)。ログ (Python
     テストが参照):
     `[silicide] nickel -> nisi 773.15 K 50 s: x 0 -> 0.0995 um
     (dSi=0.0816 um, dMet=0.0448 um, recede=0.0269 um), mesh N tets`
     (数値は `fmt("%.4g", ...)`)。
5. **デッキコマンド** (`src/deck.cpp`):
   `silicide metal=nickel time=1min temp=500C` (`metal` 必須、
   nickel/ni/titanium/ti 以外は `c.fail`)。ディスパッチ連鎖に追加。
6. **pybind11 バインディング** (`python/_cprocess.cpp`、CLAUDE.md 規約):
   ```cpp
   m.def("proc_silicide",
       [](SimState& st, const std::string& metal, double temp_k,
          double time_s) {
         std::ostringstream log;
         proc::silicide(st, metal, temp_k, time_s, &log);
         return log.str();
       },
       py::arg("state"), py::arg("metal"), py::arg("temp_k"),
       py::arg("time_s"),
       "Blanket silicidation of a deposited metal film. K/s core units.\n"
       "Returns a log string.");
   ```
7. **Simulation メソッド** (`python/cprocess/simulation.py`、oxidize_2d の
   直後):
   ```python
   def silicide(self, metal: str, time: float, temp: float) -> "Simulation":
       """Silicidation of a deposited blanket metal film (P3-d).

       metal: 'nickel' or 'titanium' (deposit it first). time in minutes,
       temp in Celsius. Diffusion-limited growth x^2 = B*t; consumes
       silicon (NiSi: 0.82, TiSi2: 0.90 of the silicide thickness) and
       metal, receding the surface. Stops when the metal is exhausted.
       """
       self._emit(_c.proc_silicide(self._st, metal,
                                   _celsius_to_k(temp), time * MIN))
       return self
   ```
8. `CMakeLists.txt` の `foreach(t ...)` に `silicide` を追加。

**ParamDB キーと既定値** (全て materials.hpp のキー一覧コメントにも記載):

| キー | 既定値 | 意味 / 根拠 |
|---|---|---|
| `silicide.nisi.b0`   | 1.2e-2 cm²/s | x(500°C, 50 s) ≈ 0.10 µm となる較正 |
| `silicide.nisi.eb`   | 1.5 eV | NiSi 拡散律速の文献オーダー |
| `silicide.nisi.rsi`  | 0.82 | 文献: Ni 1 nm → Si 1.83 nm 消費 → NiSi 2.2 nm (t_Si/t_NiSi = 1.83/2.2 ≈ 0.82) |
| `silicide.nisi.rmet` | 0.45 | 同 (t_Ni/t_NiSi = 1/2.2 ≈ 0.45) |
| `silicide.tisi2.b0`  | 5.0e-3 cm²/s | x(750°C, 600 s) ≈ 0.064 µm となる較正 |
| `silicide.tisi2.eb`  | 1.8 eV | C54 実効 1 相の文献オーダー |
| `silicide.tisi2.rsi` | 0.90 | 文献: Ti 1 nm → Si 2.27 nm 消費 → TiSi₂ 2.51 nm (2.27/2.51 ≈ 0.90) |
| `silicide.tisi2.rmet`| 0.40 | 同 (1/2.51 ≈ 0.40) |
| `<Sym>.seg_sil_m0` / `<Sym>.seg_sil_e` | 0.3 / 0.0 | m_sil = C_si/C_sil < 1: シリサイド側優先偏析 (dose loss) |
| `<Sym>.dsil0` / `<Sym>.esil` | 1.0e-3 / 2.0 | シリサイド中の実効 D (非 0 でないと偏析交換が遮断される) |

## テスト仕様

`tests/test_silicide.cpp` (新規)。基準メッシュ (薄膜分解能重視):
`mesh_box(st, 0, 0.1e-4, 0, 0.1e-4, 0, 0.5e-4, 2, 2, 250)` (z セル高 2 nm)。

1. **成長則 x ∝ √t (両対数勾配 0.5 ± 0.05)**: 独立な 3 state で
   `deposit("nickel", 0.2e-4)` (Ni 200 nm) →
   `silicide(st, "nickel", 773.15, t)`、t = {50, 200, 800} s
   (解析値 x ≈ 99.5 / 199 / 398 nm、必要金属 rmet·x ≤ 179 nm < 200 nm)。
   メッシュは z=1.0 µm, nz=500 (h=2 nm) を使用。
   (a) 戻り値の勾配 `(ln x3 − ln x1)/(ln t3 − ln t1)` が 0.5 と誤差
   < 1e-9 (解析式なので厳密)、(b) **タグから計測した**シリサイド帯厚
   (nisi セルの最大ノード z − 最小ノード z) での同勾配が **0.5 ± 0.05**
   (±1 セル ≈ ±2% の離散化を含む)。
2. **Si 消費比が文献値 ± 5%**: 基準メッシュ + Ni 60 nm →
   `silicide("nickel", 773.15, 50)` (x ≈ 99.5 nm)。silicide 前後の
   z_si_top (is_silicon 最大ノード z) の差 dSi_meas とタグ計測の
   シリサイド厚 x_meas から `dSi_meas / x_meas` が
   **0.82 ± max(0.05·0.82, 2h/x_meas)** (h=2 nm → 離散化 ±4%)。
   TiSi₂ 変種: Ti 40 nm → `silicide("titanium", 1023.15, 600)`
   (x ≈ 63.6 nm) で同様に **0.90 ± max(0.045, 2h/x_meas)**。
   さらに表面後退: bbox は不変なので「非 gas セルの最大ノード z」の低下量
   が `(rsi+rmet−1)·x` と ±2 セルで一致。
3. **金属律速キャップ**: Ni 20 nm のみ deposit → `silicide("nickel",
   773.15, 800)`。戻り値が `0.020e-4 / 0.45` cm (≈44.4 nm) と相対差
   < 1e-9、金属セルが残存しない (kMatMetal セル数 0)、layer_stack から
   金属エントリが消えている、ログに `metal fully consumed`。
4. **偏析 dose loss**: 基準メッシュ、`init("B", 1e18)` → Ni 60 nm →
   `silicide("nickel", 773.15, 50)` → `diffuse(30 min, 800 °C)`
   (DiffuseOpts, plain run で可)。
   (a) 全材料合計の B 総量 (Σ C·V) が silicide 直前比で相対差 < 2%
   (再タグのみなので保存はほぼ厳密)、(b) Si セル内 B 総量が diffuse 後に
   silicide 直後比で **> 5% 減少** (dose loss の存在)、(c) 界面両側
   1 セルの濃度比 `C_si / C_sil` が `segregation_m_silicide(B, 1073.15)`
   (= 0.3) と **factor 1.3 以内** (0.23–0.39)。
5. **異常系**: (a) 金属なし (裸 Si) で silicide → throw、(b) Ni の上に
   さらに oxide を deposit した state で throw、(c) 不明金属 "cu" →
   throw、(d) `set_param("silicide.nisi.rsi", ...)` 等のキャリブレーション
   が効く (rsi=0.5 に上書きして dSi/x ≈ 0.5 ± 離散化)。
6. **デッキ**: `run_deck` に `"silicide metal=nickel time=1min temp=500C"`
   を含む小デッキ (deposit nickel 込み) が完走し `[silicide]` を含む。

Python (`python/test_comprehensive.py` に `test_silicide_python()` を追加):

1. `sim.mesh(0.1, 0.1, 0.5, 2, 2, 250).init("B", 1e18)
   .deposit("nickel", 0.06).silicide("nickel", time=50/60, temp=500)` が
   通り、(a) チェーン (`self` 返し) 成立、(b) ログに `[silicide]`、
   (c) 非 gas 最上面が deposit 直後より低い (表面後退、cell_centroids で
   確認)、(d) `sim.diffuse(30, 800)` 後、Si 帯の B 総量が減少 (numpy)。
2. `deposit("titanium", 0.04).silicide("titanium", 10, 750)` の smoke
   (完走 + シリサイド厚 > 0 をログ数値で確認)。
3. 既存の Python テストが全て PASS。

既存 C++ テスト (`ctest --test-dir build`) が全て PASS すること。
SegTable/MatId 拡張は共有コードなので特に `test_segregation`,
`test_multimaterial`, `test_oed`, `test_oxidize_flow`, `test_ox2d`,
`test_diffusion`, `test_state_io` (region_material 文字列の保存/復元に
"nickel"/"nisi" が乗るケースも 1 ケース追加) の回帰を必ず確認。

## 完了条件 (DoD)

- [ ] MatId に kMatSilicide/kMatMetal 追加、SegTable が kMatCount=7 次元、
      material_id / material_diffusivity / make_seg_table 拡張
- [ ] Dopant に seg_sil_m0/seg_sil_e/dsil0/esil カラム (末尾追加) +
      アクセサ + ParamDB キー (上表の 12 キー) + ヘッダのキー一覧更新
- [ ] `proc::silicide` 実装 (再タグのみ、x²=x₀²+B·t、金属律速キャップ、
      layer_stack 維持、戻り値 = 総シリサイド厚 cm)
- [ ] deposit/set_region が nickel/titanium を受理
- [ ] デッキコマンド `silicide` (metal/time/temp) 追加
- [ ] `proc_silicide` バインディング + `Simulation.silicide` (min/°C 変換、
      `self` 返し) — CLAUDE.md の必須 3 点セット
- [ ] `tests/test_silicide.cpp` の 6 テスト PASS、`CMakeLists.txt` の
      foreach に `silicide` 追加
- [ ] `test_comprehensive.py` に Python テスト追加、全 Python テスト PASS
- [ ] `cmake --build build -j$(nproc) && ctest --test-dir build` 全 PASS
- [ ] コミットメッセージに `P3-d` を含める

## やらないこと

- **TiSi₂ の C49→C54 相転移カイネティクス** — 実効 1 相 (単一の B(T))。
  Ni₂Si→NiSi→NiSi₂ の相列も同様に NiSi 実効 1 相のみ
- 核生成律速・線形 (反応律速) レジーム — 純拡散律速 x²=B·t のみ
  (Deal-Grove の線形項 A は導入しない; deal_grove_step は流用しない)
- **横方向侵食 (lateral encroachment)・2D 形状** — ブランケット 1D のみ。
  ゲート脇のスペーサ選択性等は将来 oxidize_2d 型の別エントリポイント
- ポリサイド (poly 上のシリサイド化) — 表面は単結晶 Si のみ許容 (throw)
- 金属/シリサイドの MC 注入ターゲット材追加 (P3-a の成分 BCA の範囲)
- 体積変化に伴う応力 (P2-6/P3-e)、Kirkendall ボイド、シート抵抗出力
- スノープラウ (成長フロントでの非平衡押し出し) — dose loss は
  SegTable の平衡偏析 + 界面輸送 h のみ
- メッシュの物理的縮小 — 表面後退は gas 再タグで表現 (polygon etch と
  同じ制約; 真のセル除去は P2-5 系の一般化スコープ)
