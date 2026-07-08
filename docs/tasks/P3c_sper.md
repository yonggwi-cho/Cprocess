# P3-c: 固相エピ再成長 (SPER)

## 目的

MC 注入でアモルファス化した表面領域 (変位原子密度 damage ≥ 6.25e21 cm⁻³ =
`kAmorphizationDensity`) は、低温アニール (550–650 °C) で結晶/アモルファス
界面 (a/c フロント) が下から上へ掃引して c-Si に再成長する。本タスクで:

- アモルファス深さ z_a(x,y) を **damage 場の閾値交差**からカラム毎に抽出し、
  フロント速度 v = v0·exp(−Ea/kT) (v0 = 3.1e8 cm/s、Ea = 3.1 eV、<100>) で
  z_a をカラム毎 1D 時間発展させる
- 再成長した領域の damage を 0 化し (TED の I シードも同時に消す)、
  **SPER 準安定活性化** (固溶度 C_ss(T) の up to `sper.act_factor` 倍、
  既定 10) を与える
- `proc::sper(st, temp_k, time_s, log)` を **diffuse_ted の前処理として使う
  独立コマンド**にする (P3 概要の確定事項。スノープラウは第一段階で省略)

**P1-2 (MC 損傷 → I シード) と P2-1 (フル点欠陥) はマージ済みで前提。**

本仕様の確定定数から導かれる代表値 (テストの数値根拠、kB = 8.617333262e-5
eV/K で計算):

| T | kT [eV] | v = 3.1e8·exp(−3.1/kT) [cm/s] | 100 nm 全再成長時間 |
|---|---------|-------------------------------|---------------------|
| 550 °C (823.15 K) | 0.0709335 | 3.25e-11 | 3.08e5 s |
| 600 °C (873.15 K) | 0.0752422 | 3.97e-10 | **2.52e4 s (≈ 420 min)** |
| 650 °C (923.15 K) | 0.0795510 | 3.69e-9  | 2.71e3 s (≈ 45 min) |

概要の「600 °C で ~分オーダー」は確定定数 (Ea=3.1 eV, v0=3.1e8) からは
10² 分オーダー (420 min) になる。**定数が確定事項なので、受入テストは
計算値 2.52e4 s を採用する** (650 °C なら文字どおり分オーダー ≈ 45 min)。

## 現状コード

- `src/mc_implant.cpp` — `apply_mc_implant(..., damage_conc)`: channeling 有効時
  のみ Kinchin-Pease 変位原子密度 [cm⁻³] を `damage_conc` に返す (L607-612)。
  `DamageStore::f_amor` は `dens / kAmorphizationDensity` で結晶を劣化させる。
- `src/process.cpp` — `implant_mc` (L302): 生の damage 配列 `dmg` を
  `seed_interstitials_from_damage(st, dmg)` (L137) に渡して
  `st.fields["I"] += kFrenkelSurvival * damage` (上限 `kAmorphizationDensity`)
  したあと**破棄する**。つまり P1-2 の「損傷場」は 0.01 倍・キャップ済みの
  "I" 場としてしか永続化されておらず、**生の変位原子密度そのものは
  st.fields に残らない** → 本タスクで `st.fields["damage"]` として永続化する
  (下記手順 1)。stack 経由ブランチ (L348-351) と直接ブランチ (L380-381) の
  2 箇所がシード地点。
- `include/cprocess/materials.hpp` — `kFrenkelSurvival = 0.01` (L49)、
  `kAmorphizationDensity = 6.25e21` (L55)、`solid_solubility(d, T)`
  (materials.cpp L173、ParamDB `<Sym>.ss_pre / <Sym>.ss_e`)、
  `active_concentration(d, conc, T) = min(C, C_ss)`。
  As: ss_pre=1.3e23, ss_e=0.66 → **C_ss(As, 600 °C) = 2.02e19 cm⁻³**。
- `src/process.cpp` — `proc::active_field` (L2150) と `proc::save` (L2115、
  `<Sym>_active` / `NetDoping` の組立) が `active_concentration` を呼ぶ 2 つの
  活性化クランプ地点。SPER 準安定活性化はこの 2 箇所に反映する (手順 4)。
- `src/process.cpp` — `proc::diffuse_ted` (L1862): `st.fields["I"/"V"/"C311"]`
  を作成し `DiffusionSolver::run_ted(fields, psi, v, c311, bcs, opts)` (P2-1
  の 4 場版) を呼ぶ。sper は run_ted に**手を入れない** (前処理コマンド)。
- `include/cprocess/param_db.hpp` — `ParamDB::instance().get(key, fallback)`
  シングルトン。materials.cpp と同じパターンで全定数を読む。
- `src/state_io.cpp` — CPRC1 は st.fields を名前つきで全部シリアライズする
  ので、新フィールド "damage" / "regrown" の save/load 対応は**自動**。
- `src/process.cpp` — `extend_mesh_exact` / `infer_box_dims` (L101-108):
  box メッシュの (nx, ny, nz) を復元できる。カラム分解はこれを使う。
- `python/_cprocess.cpp` / `python/cprocess/simulation.py` — proc_* バインディング
  と `Simulation` メソッドの定型 (CLAUDE.md の必須ラッパ規約)。

## 実装手順

1. **生 damage 場の永続化 (process.cpp)**。無名 namespace に追加:
   ```cpp
   // Persist the raw MC displaced-atom density [cm^-3] (uncapped, no
   // kFrenkelSurvival factor) so SPER (P3-c) can find the amorphous depth.
   void accumulate_damage(SimState& st, const std::vector<double>& dmg) {
     auto& D = st.fields["damage"];
     D.resize(st.mesh.cells.size(), 0.0);
     for (std::size_t i = 0; i < dmg.size() && i < D.size(); ++i) D[i] += dmg[i];
   }
   ```
   `implant_mc` の `use_damage` 分岐 2 箇所 (stack: L351 の
   `seed_interstitials_from_damage(st, dmg_transferred)` 直後に
   `accumulate_damage(st, dmg_transferred)`、直接: L381 直後に
   `accumulate_damage(st, dmg)`) に挿入。"I" のシード則・キャップは不変。
   `implant_gauss` / "+1" モデルでは damage 場を作らない (MC 専用)。

2. **ParamDB キー** (materials.hpp のキー表コメントに追記):

   | キー | 意味 | 既定値 |
   |------|------|--------|
   | `sper.v0` | フロント速度前指数 [cm/s] (<100>) | 3.1e8 |
   | `sper.ea` | 活性化エネルギー [eV] | 3.1 |
   | `sper.amorph_density` | アモルファス化閾値 [cm⁻³] | `kAmorphizationDensity` (6.25e21) |
   | `sper.act_factor` | 準安定活性化の C_ss 倍率 | 10.0 |

   すべて呼び出し時に `ParamDB::instance().get(key, default)` で読む。

3. **`proc::sper` (process.hpp / process.cpp)**:
   ```cpp
   // P3-c: solid-phase epitaxial regrowth. Isothermal anneal at temp_k for
   // time_s. Reads st.fields["damage"] (raw MC displaced-atom density,
   // persisted by implant_mc); columns whose surface-connected damage exceeds
   // sper.amorph_density regrow bottom-up at v = sper.v0*exp(-sper.ea/kT).
   // Regrown cells: damage = 0, I = 0, fields["regrown"] = 1 (metastable
   // activation marker consumed by active_field()/save()). Intended to run
   // BEFORE diffuse_ted. No-op (with log) if no damage field / no amorphous
   // cells. Throws if a photoresist stack is present.
   void sper(SimState& st, double temp_k, double time_s,
             std::ostream* log = nullptr);
   ```
   実装 (すべて Si セルのみ、`silicon_mask(st)`):
   1. `need_mesh(st)`; `st.has_stack` なら throw ("strip before sper");
      `time_s < 0 || temp_k <= 0` は throw。`st.fields` に "damage" が無い、
      または `max(damage) < c_am` (`c_am = get("sper.amorph_density", ...)`)
      なら `[sper] no amorphous region (nothing to do)` をログして return。
   2. **カラム分解**: `infer_box_dims(st.mesh)` で (nx, ny) を得て
      `hx = (bb.hi.x-bb.lo.x)/nx`, `hy` 同様。セル ci の xy-bbox 最小値から
      `ix = llround((xmin - bb.lo.x)/hx)`, `iy` を計算し `key = iy*nx + ix` で
      グループ化 (box メッシュでは 1 ブリックの 5-6 tet が同じ key を持つ。
      profile1d と同趣旨の xy-bbox 判定の格子版)。
   3. カラム毎に Si セルを centroid.z **降順** (表面→深部) にソート。
      `z_top` = カラム内 Si セル bbox の最大 z。表面セルから damage ≥ c_am
      が連続する区間 (surface-connected run) を取り、その最深セルの
      bbox 最小 z を `z_bot0` として初期アモルファス厚
      `d_a0 = z_top − z_bot0` とする。表面セルが非アモルファスのカラムと、
      連続区間より深い埋め込みアモルファス帯はスキップ (件数をログ)。
   4. **時間発展 (解析)**: 等温なので v は定数。
      `v = get("sper.v0",3.1e8) * exp(−get("sper.ea",3.1)/(kBoltzmannEv*temp_k))`、
      再成長厚 `dz = min(v * time_s, d_a0)`、残アモルファス厚
      `d_a = d_a0 − dz` (フロントは a/c 界面 z_bot0 から上向きに進む)。
      サブステップは不要 (温度ランプは対象外、やらないこと参照)。
   5. **再成長セルの判定**: 初期アモルファス区間内のセルで
      `centroid.z < z_top − d_a` (= フロントが centroid を通過) のものを
      再成長とする。再成長セルで:
      - `st.fields["damage"][ci] = 0.0`
      - `st.fields["I"]` が存在すれば `I[ci] = 0.0` (アモルファス内の
        カスケード損傷はフロント通過で消滅。EOR 側 = 区間より深いセルの
        I は不変 → その後の diffuse_ted の TED 源として残る)
      - `st.fields["regrown"]`(無ければ n_cells で 0 初期化) `= 1.0`
   6. ログ (1 行 + 集計):
      `[sper] T=<temp> K, t=<time> s, v=<v> cm/s (Ea=<ea> eV)` と
      `[sper] columns=<n>, z_a max <d_a0max*1e4> -> <d_amax*1e4> um, regrown cells=<k>` 。
   7. `st.last_temp` は**更新しない** (拡散を伴わない。活性化温度は
      active_field 側で明示指定させる)。

4. **準安定活性化の反映 (process.cpp のみ、materials.cpp は不変)**:
   `proc::active_field` (L2150) と `proc::save` の active/NetDoping ループ
   (L2124-2134) の 2 箇所で、`active_concentration(*d, conc[i], t_k)` を
   次に置き換える:
   ```cpp
   const double css = solid_solubility(*d, t_k);
   const double f_sper = ParamDB::instance().get("sper.act_factor", 10.0);
   const auto* rg = /* st.fields.find("regrown") の結果、無ければ nullptr */;
   const double cap = (css > 0 && rg && (*rg)[i] > 0.5) ? css * f_sper : css;
   act[i] = (cap > 0) ? std::min(conc[i], cap) : conc[i];
   ```
   (css/f_sper/rg はループ外で 1 回だけ引く。) "regrown" 場は state_io で
   自動永続化されるので save_state/load_state 越しに活性化が保たれる。

5. **pybind11 バインディング** (python/_cprocess.cpp、proc_* 規約どおり):
   ```cpp
   m.def("proc_sper",
       [](SimState& st, double temp_k, double time_s) {
         std::ostringstream log;
         proc::sper(st, temp_k, time_s, &log);
         return log.str();
       },
       py::arg("state"), py::arg("temp_k"), py::arg("time_s"));
   ```

6. **`Simulation.sper`** (python/cprocess/simulation.py、diffuse/oxidize の
   隣に配置):
   ```python
   def sper(self, temp: float, time: float) -> "Simulation":
       """Solid-phase epitaxial regrowth (P3-c). temp in Celsius, time in
       minutes. Regrows the MC-amorphized layer (fields 'damage' >=
       sper.amorph_density) column by column at v = v0*exp(-Ea/kT); regrown
       cells get damage/I zeroed and metastable activation (sper.act_factor x
       C_ss). Run before diffuse(ted=True)."""
       self._emit(_c.proc_sper(self._st, _celsius_to_k(temp), time * MIN))
       return self
   ```
   単位変換 (°C→K, min→s) と委譲のみ。プロセスロジックは書かない。

7. **テスト** (次節): `tests/test_sper.cpp` 新規、CMakeLists.txt L76 の
   `foreach(t ... rcm)` に `sper` を追加。Python テストを
   `python/test_comprehensive.py` に追加。

## テスト仕様

`tests/test_sper.cpp` (test_util.hpp の CHECK を使用):

1. **Arrhenius フィット誤差 < 5%**: mesh_box 0.1×0.1×1.0 µm (2×2×200、
   セル高 5 nm)。`st.fields["damage"]` を直接セット: 深さ ≤ 0.5 µm の
   セルに 1e22、他 0。T ∈ {823.15, 873.15, 923.15} K それぞれ独立の state
   で `proc::sper(st, T, t_T)`、t_T は解析 v で再成長 300 nm になる値
   {9.23e5, 7.56e4, 8.13e3} s。実測 v_meas(T) = (0.5 µm − 残アモルファス厚)
   / t_T (残厚 = damage ≥ 6.25e21 のカラム内連続セル数 × 5 nm)。
   - 各点で `|v_meas − v_analytic| / v_analytic < 0.02` (セル量子化
     5 nm/300 nm ≈ 1.7% を許容)
   - 3 点最小二乗で `ln v = ln v0 − Ea/(kB T)` をフィットし
     **`|Ea_fit − 3.1| / 3.1 < 0.05`** を CHECK
2. **100 nm 全再成長の時間窓 (600 °C)**: 同メッシュ、深さ ≤ 0.1 µm に
   damage = 1e22。解析全再成長時間 t_full = 1e-5 / 3.966e-10 = **2.522e4 s**。
   - `sper(873.15 K, 2.60e4 s)` → 全セルで damage < 6.25e21 (完全再成長)
   - 新しい同一 state で `sper(873.15 K, 2.40e4 s)` → damage ≥ 6.25e21 の
     セルが残る (未完)。残厚が `0.1 µm − v·t = 47.7 Å` オーダー…ではなく
     `1e-5 − 3.966e-10*2.40e4 = 4.8e-7 cm = 4.8 nm` ± 1 セル (5 nm) 以内
3. **As 準安定活性化 > C_ss(600 °C)**: 上記 100 nm アモルファス state に
   `st.fields["As"] = 1e20` (全セル)。`sper(873.15, 2.60e4)` 後
   `proc::active_field(st, "As", 873.15)`:
   - 再成長セル (regrown=1): act = 1e20 (**> C_ss(As,600 °C) = 2.02e19**、
     1e20 < 10×2.02e19 なのでクランプされない。相対差 < 1e-12)
   - 非再成長の深部セル: act = 2.02e19 ± 1% (通常クランプ)
   - `proc::set_param(st, "sper.act_factor", 1.0)` 後は再成長セルも
     act = 2.02e19 ± 1%
4. **I の 0 化と TED 前処理**: 手順 2 の state で "I" に 1e20 をセットして
   から全再成長 → 再成長セルで I == 0 (bit)、アモルファス区間より深い
   セルの I は不変 (bit)。続けて B=1e18 を置いて
   `proc::diffuse_ted(st, {time=60, temp=1173.15})` が例外なく完走。
5. **no-op / エラー系**: damage 場なしの state で `sper(873.15, 600)` →
   throw せず全 fields bit 不変、ログに "no amorphous"。`time_s < 0` は
   throw。`photo` 後 (has_stack) は throw、`strip` 後は成功。
6. **MC 統合**: mesh_box 0.2×0.2×0.5 µm (4×4×50) → `implant_mc("As",
   dose=1e15, 30 keV, ions=20000, channeling=true, seed_damage=true)` →
   `st.fields["damage"]` が存在し max > 0。max ≥ 6.25e21 の場合、
   `sper(873.15, 3.0e4)` でアモルファスセル数が減少する。

`python/test_comprehensive.py` に追加:

- `sim.mesh(...).init("B", 1e15)` → `implant(species="As", ..., mc=True,
  channeling=True, damage=True)` → `"damage" in sim.field_names()` →
  `ret = sim.sper(600, 500)` が完走し `ret is sim` (チェーン) →
  `sim.field("damage")` の max がsper 前より減少または不変 →
  `sim.set_param("sper.act_factor", 1.0)` 後 `sim.active("As", 600)` の
  max が `solid_solubility` 相当 (≤ 2.1e19) に戻ること (アモルファス化した
  場合のみの条件付き assert で可)。

**既存テストの回帰**: `test_mc`, `test_ted`, `test_activation`,
`test_dopant_clusters`, `test_state_io`, `test_integration`, `test_flow` が
無変更で PASS すること (implant_mc への追記は damage 場の追加のみで
"I" シード・dose は不変; active の式変形は "regrown" 場が無い限り
`active_concentration` と同値)。`python3 python/test_comprehensive.py` PASS。

## 完了条件 (DoD)

- [ ] `implant_mc` が生 damage を `st.fields["damage"]` に永続化 (両ブランチ)、
      state_io 経由で save/load される
- [ ] `proc::sper` 実装 (カラム毎 z_a 抽出 + 解析時間発展 + damage/I 0 化 +
      "regrown" マーキング)、全定数が ParamDB キー (`sper.v0` / `sper.ea` /
      `sper.amorph_density` / `sper.act_factor`) 経由
- [ ] `active_field` / `save` が regrown セルで C_ss × act_factor を上限に
      使う (準安定活性化)
- [ ] pybind `proc_sper` + `Simulation.sper(temp, time)` (°C/min 変換、
      `self` 返し) + Python テスト (CLAUDE.md 必須ラッパ規約の 3 点セット)
- [ ] `tests/test_sper.cpp` 6 テスト PASS、CMakeLists の foreach に `sper`
      登録、既存テスト全 PASS (`ctest --test-dir build`)
- [ ] `build/_cprocess*.so` を `python/cprocess/` にコピーし
      `python3 python/test_comprehensive.py` PASS
- [ ] コミットメッセージに `P3-c` を含める

## やらないこと

- **スノープラウ** (フロントでの不純物押し出し・偏析) — 概要の確定事項、
  第一段階では省略 (フロント通過で濃度場は不変)
- 温度ランプ中の SPER (等温のみ。RTA ランプ対応は将来 DiffuseOpts 連携で)
- 面方位依存 (<110>/<111> の v0 縮小)、ドーパント濃度・不純物 (O, F) による
  v の増減速 — v は温度のみの関数
- 埋め込み (表面非連結) アモルファス帯の再成長 (両面フロント)。表面連結の
  単一区間のみ、フロントは 1 本
- ランダム核形成 (polycrystallization)、EOR 転位ループの明示モデル
  (既存 "I"/"C311" 場に委ねる)
- 準安定活性化のその後の高温アニールでの脱活性化 (regrown フラグは
  立ちっぱなし。deactivation kinetics は将来タスク)
- run_ted / diffusion.cpp 内部の変更 (活性化ゲート diffusion.cpp:1128 への
  sper_cap 反映を含む — 反映先は active_field/save のみ)
- テキストデッキコマンド (`sper temp=... time=...`) の追加 (Python 経由で
  十分。必要なら別タスク)
- `implant_gauss` (解析注入) からの damage 生成 — MC + channeling 専用
