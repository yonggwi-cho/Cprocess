# P3-b: エピタキシー成長 (`proc::epitaxy`)

## 目的

その場ドープを含む Si エピタキシャル層の成長を 1 コマンド化する。
幾何は deposit の機構 (extend_mesh_exact + 再タグ + layer_stack) を材料
"silicon" で流用し、新規セルには doping_map (種→cm⁻³) の一様濃度を設定、
成長温度での熱履歴として成長時間分の diffuse を自動で 1 回実行する
(`anneal=true` 既定)。成長速度モデル (H₂/SiH₄ 律速) は導入せず、
thickness/time は直接指定とする (P3_overview.md の確定方針)。

**P1-7 (etch/depo + layer_stack) と M-2 (adaptive refine) の完了が前提
(両方マージ済み)。**

## 現状コード

- `src/process.cpp` — `deposit(st, material, thickness, nz_add, poly, log)`
  (L502): 材料 "silicon"/"si" を既にサポート (`known[]` リスト L510)。
  `extend_mesh_exact(st.mesh, thickness, nz_add)` で上方拡張し、
  合成タグ `dep_tag = max_tag + 1000` を新設、centroid > z_top の新規セルを
  `dep_tag` に、既存セルは最近傍 centroid で region/全フィールドを転写
  (**新規セルのドーパントは 0 に初期化される** — 本タスクはこの直後に
  doping_map で上書きする)。最後に
  `st.layer_stack.insert(st.layer_stack.begin(), {dep_tag, mat})`。
  **epitaxy はこの deposit() をそのまま内部呼び出しして幾何を実現する**
  (再実装しない)。
- `src/process.cpp` — 無名 namespace の `extend_mesh_exact` (L101):
  `infer_box_dims` で nx/ny/nz を復元し `make_box_mesh` で作り直す。
  bbox 上面は厳密に `+thickness` されるので、エピ層厚は離散化誤差なしで
  実現される (z 分割は全体で均等再分割、よって既存セルは最近傍転写)。
- `src/process.cpp` — `oxidize()` の OED パス (L1480–1493) が
  「プロセスステップ内から `proc::diffuse_ted` を自動実行する」手本:
  `DiffuseOpts opts_k; opts_k.temp = temp_k; opts_k.time = dt_k;
  opts_k.verbosity = 0; opts_k.lin_maxit = 5000; opts_k.lin_rtol = 1e-8;`
  (新鮮な材料界面直後の連成系ではデフォルトの lin_rtol=1e-10/2000 が
  収束しないことがあるため緩和済み — epitaxy も同じ設定を使う)。
- `src/process.cpp` — `diffuse_ted(st, opts, log)` (L1862): ドーパント場が
  無ければ `[ted] no dopants present` を出して安全に return。"I" 場が
  無い/ゼロなら平衡アニール (S=1) に退化する — つまり **implant 履歴の
  有無によらず常に diffuse_ted を呼んでよい**。
- `src/process.cpp` — `init()` (L227): `silicon_mask` でマスクした一様濃度
  設定の手本 (ただし本タスクは dep_tag セル限定なので直接ループする)。
- `include/cprocess/deck.hpp` — `SimState::layer_stack`
  (`std::vector<std::pair<int,std::string>>`、最新層が先頭)。
- `src/deck.cpp` — `cmd_oxidize` が time/temp の `Unit::` 換算 +
  ディスパッチ追加の手本。
- Python: `python/_cprocess.cpp` の `proc_etch_rate` が
  `std::map<std::string,double>` 引数を受けるバインディングの手本 (L~560)、
  `python/cprocess/simulation.py` の `deposit()` (L214) / `oxidize()`
  (L292) がメソッドの手本 (`UM = 1e-4`, `MIN = 60.0`, `_celsius_to_k`)。
- ParamDB (P1-10): **本タスクは新しい物理定数を導入しない** (厚さ/時間
  直接指定 + 既存拡散パラメータのみ)。新規 ParamDB キーはなし。

## 実装手順

1. `include/cprocess/process.hpp` に宣言を追加 (deposit_conformal の直後):
   ```cpp
   // Epitaxial growth of `thickness_cm` of silicon on the exposed silicon
   // top surface, with in-situ uniform doping `doping` (species symbol ->
   // cm^-3) in the newly grown cells only. Geometry reuses deposit()'s
   // extend_mesh_exact + retag + layer_stack mechanism with material
   // "silicon". When anneal == true (default) the growth thermal budget is
   // applied by one automatic proc::diffuse_ted(temp_k, time_s) call after
   // growth, so substrate dopants back-diffuse into the epi layer. No
   // growth-rate model: thickness and time are both caller-given.
   void epitaxy(SimState& st, double thickness_cm, double temp_k,
                double time_s,
                const std::map<std::string, double>& doping = {},
                bool anneal = true, std::ostream* log = nullptr);
   ```
2. `src/process.cpp` に実装 (すべて cm/s/K):
   - (a) `need_mesh(st)`; `thickness_cm <= 0` / `time_s <= 0` は throw;
     `st.has_stack` なら `throw std::runtime_error("epitaxy: strip resist
     first")`。doping の各キーを `dopant_or_throw(sym)` で検証 (不明種は
     ここで throw、幾何変更前に fail-fast)。
   - (b) **表面が Si であることの検証** (oxidize L1181–1198 の走査と同型):
     `z_top = st.mesh.bbox().hi.z`、z セル高
     `h = (bbox z 幅) / nz` (`infer_box_dims`)。centroid が
     `z > z_top - h` の全セルについて材料 (`region_material`、未タグは
     silicon) が `is_silicon` でなければ
     `throw std::runtime_error("epitaxy: top surface is not silicon")`
     (酸化膜/窒化膜上の poly 核生成はモデル外 — 「やらないこと」参照)。
   - (c) **成長**: `nz_add = std::max(1, (int)std::round(thickness_cm / h))`
     として `deposit(st, "silicon", thickness_cm, nz_add, {}, nullptr);`
     を呼ぶ (ログは epitaxy 自身が 1 行で出すので deposit のログは抑制)。
     エピ層タグは `const int epi_tag = st.layer_stack.front().first;`
     (deposit が直前に push_front した項目)。
   - (d) **その場ドープ**: 各 `(sym, conc)` について
     ```cpp
     auto& f = st.fields[find_dopant(sym)->symbol];   // 正規化シンボルで格納
     f.resize(st.mesh.cells.size(), 0.0);
     for (std::size_t i = 0; i < f.size(); ++i)
       if (st.mesh.cell_region[i] == epi_tag) f[i] = conc;  // 一様、上書き
     ```
     conc < 0 は (a) で throw。conc == 0 は許容 (明示アンドープ)。
   - (e) **成長熱履歴のアニール**: `anneal == true` のとき
     ```cpp
     DiffuseOpts opts;
     opts.temp = temp_k;  opts.time = time_s;  opts.verbosity = 0;
     opts.lin_maxit = 5000;  opts.lin_rtol = 1e-8;   // oxidize(P2-3) と同じ緩和
     diffuse_ted(st, opts, log);
     ```
     を **1 回だけ** 呼ぶ (P3_overview の確定方針: サブステップ分割なし)。
     diffuse_ted は "I" 場が無ければ平衡アニールに退化し、ドーパント場が
     無ければ no-op なのでガード不要。`anneal == false` は幾何+ドープのみ。
   - (f) `st.last_temp = temp_k;` (anneal の有無によらず — 成長自体が
     temp_k での熱工程)。
   - (g) ログ (1 行、Python テストが `[epitaxy]` を参照):
     `[epitaxy] 0.1 um Si @ 1273.15 K 600 s, doping: P=1e+18, anneal=on,
     mesh now N tets` (doping 空なら `undoped`、値は `fmt("%.4g", ...)`)。
3. **デッキコマンド** (`src/deck.cpp`、`cmd_oxidize` に倣う):
   ```
   epitaxy thickness=0.5um temp=1000C time=10min [species=B conc=1e17] [anneal=on|off]
   ```
   `thickness` は `Unit::len`、`temp`/`time` は既存換算。デッキでは
   単一種ドープのみ (species/conc ペア、両方あるか両方ないか;
   片方だけなら `c.fail(...)`)。多種ドープは Python 専用と README 級の
   コメントで明記。`run_deck` のディスパッチ連鎖に追加。
4. **pybind11 バインディング** (`python/_cprocess.cpp`、CLAUDE.md 規約):
   ```cpp
   m.def("proc_epitaxy",
       [](SimState& st, double thickness_cm, double temp_k, double time_s,
          const std::map<std::string, double>& doping, bool anneal) {
         std::ostringstream log;
         proc::epitaxy(st, thickness_cm, temp_k, time_s, doping, anneal, &log);
         return log.str();
       },
       py::arg("state"), py::arg("thickness_cm"), py::arg("temp_k"),
       py::arg("time_s"), py::arg("doping") = std::map<std::string, double>{},
       py::arg("anneal") = true,
       "Epitaxial Si growth with in-situ doping. cm/K/s core units.\n"
       "Returns a log string.");
   ```
5. **Simulation メソッド** (`python/cprocess/simulation.py`、oxidize() の
   直前、deposit 系の並び):
   ```python
   def epitaxy(self, thickness: float, temp: float, time: float, *,
               doping: dict | None = None, anneal: bool = True) -> "Simulation":
       """Epitaxial silicon growth with in-situ doping.

       thickness in micrometres, temp in Celsius, time in minutes.
       doping: {species: conc_cm3} set uniformly in the new layer only.
       anneal=True (default) applies the growth thermal budget as one
       automatic TED-capable anneal (diffuse_ted) after growth.
       """
       self._emit(_c.proc_epitaxy(self._st, thickness * UM,
                                  _celsius_to_k(temp), time * MIN,
                                  dict(doping or {}), bool(anneal)))
       return self
   ```
6. `CMakeLists.txt` の `foreach(t ...)` に `epitaxy` を追加。

## テスト仕様

`tests/test_epitaxy.cpp` (新規)。共通メッシュ:
`mesh_box(st, 0, 0.2e-4, 0, 0.2e-4, 0, 0.5e-4, 4, 4, 50)` (z セル高 0.01 µm)。

1. **厚さ・ドーピングの一致** (`anneal=false` で幾何+ドープを単離):
   `epitaxy(st, 0.1e-4, 1273.15, 600, {{"P", 1e18}}, false)` の後、
   (a) `bbox().hi.z` の増分が 0.1 µm と相対差 < 1e-9
   (extend_mesh_exact は厳密なので **< 1 セル高どころか機械精度**)、
   (b) エピタグ (`layer_stack.front()`) の全セルで P == 1e18
   (相対差 < 1e-12、受入基準の < 1% を厳密側で満たす)、
   (c) 旧表面 (0.5 µm) 以下のセルの P は全て 0。
2. **基板ドーパントの遡り拡散が 2√(Dt) と整合**: 基板に
   `init(st, "B", 1e19)` (1373 K で ni≈1.2e19 なのでほぼ真性、電界増速の
   影響が小さい)。`epitaxy(st, 0.2e-4, 1373.15, 7200, {}, true)`
   (1100 °C 120 min、√(Dt) = √(D_B(1373K)·7200) ≈ 73 nm)。
   界面 (z=0.5 µm) からの step 拡散は C(z) = (C0/2)·erfc(z/2√(Dt)) なので
   C = C0/4 となる高さは z₁/₄ = 0.954·√(Dt) ≈ 70 nm。テスト:
   `profile1d(st, "B", 0.1e-4, 0.1e-4)` から z > 0.5 µm で最初に
   B < 0.25e19 となるセル centroid 高さ z_meas を取り、
   `|z_meas - 0.954*sqrt(D*t)| <= max(0.30 * 0.954*sqrt(D*t), 2*h)`
   (D は `dopant_diffusivity(*find_dopant("B"), 1373.15, 1.0)` で
   テスト内計算、h = 0.01 µm)。加えてエピ層内の B 総量が > 0
   (遡り拡散の存在)。
3. **多層エピ (3 回連続) の layer_stack**: 新しい state で
   `epitaxy(0.05e-4, T, t, {{"B",1e16}}, false)` →
   `{{"B",1e17}}` → `{{"B",1e18}}` の 3 連続。
   (a) `layer_stack.size()` が 3 増え、先頭 3 項目の材料が全て
   "silicon" で region タグが互いに異なる、(b) `bbox().hi.z` 増分が
   0.15 µm (相対差 < 1e-9)、(c) `profile1d` の B が上から
   1e18 / 1e17 / 1e16 の 3 段 (各層中央セルで相対差 < 1e-12)。
4. **異常系**: (a) `deposit(st, "oxide", 0.02e-4)` 後の epitaxy →
   `std::runtime_error` ("not silicon")、(b) doping に不明種
   `{"Xx", 1e18}` → throw (幾何は不変: セル数が前後で一致)、
   (c) `thickness <= 0` / `time <= 0` → throw。
5. **デッキ**: `run_deck` に
   `"epitaxy thickness=0.1um temp=1000C time=10min species=P conc=1e18"`
   を含む小デッキを食わせ、完走しログに `[epitaxy]` が含まれること。

Python (`python/test_comprehensive.py` に `test_epitaxy_python()` を追加):

1. `sim.mesh(0.2, 0.2, 0.5, 4, 4, 50).init("B", 1e15)
   .epitaxy(0.1, 1000, 5, doping={"P": 1e18}, anneal=False)` が通り、
   (a) メソッドチェーン (`self` 返し) が成立、
   (b) `sim.cell_centroids[:, 2].max()` が 0.5 µm 超、
   (c) `sim.field("P")[z > 0.5 µm のセル]` の最小値 == 1e18 (rel < 1e-9)、
   (d) ログに `[epitaxy]` を含む。
2. 同一フローの `anneal=True` 版で、界面をまたぐ B プロファイルの
   エピ層側総量が anneal=False 版より大きい (遡り拡散の存在、numpy 比較)。
3. 既存の Python テストが全て PASS。

既存 C++ テスト (`ctest --test-dir build`) が全て PASS すること。
特に `test_etch_depo` (deposit/layer_stack の回帰) と `test_ted` /
`test_oed` (diffuse_ted 呼び出し規約) は本タスクの流用元なので必ず確認。

## 完了条件 (DoD)

- [ ] `proc::epitaxy` 実装 (process.hpp/process.cpp、deposit() 内部呼び出し
      + epi_tag 一様ドープ + anneal 時の diffuse_ted 1 回)
- [ ] 新規 ParamDB キーなし (新定数を導入しないことが仕様)
- [ ] デッキコマンド `epitaxy` (thickness/temp/time/species/conc/anneal) 追加
- [ ] `proc_epitaxy` バインディング + `Simulation.epitaxy` (µm/min/°C 変換、
      `self` 返し) — CLAUDE.md の必須 3 点セット
- [ ] `tests/test_epitaxy.cpp` の 5 テスト PASS、`CMakeLists.txt` の
      foreach に `epitaxy` 追加
- [ ] `test_comprehensive.py` に Python テスト追加、全 Python テスト PASS
- [ ] `cmake --build build -j$(nproc) && ctest --test-dir build` 全 PASS
- [ ] コミットメッセージに `P3-b` を含める

## やらないこと

- **成長速度モデル** (H₂/SiH₄ 律速、温度/圧力/分圧依存レート) —
  thickness/time 直接指定 (P3_overview の確定方針)
- 選択エピ (マスク開口内のみの成長)・ファセット形成 — ブランケットのみ
  (deposit(poly=) を意図的に露出しない; 非 Si 表面は throw)
- SiGe / Ge 組成場・歪み (P3-f)。doping に "Ge" を渡すこと自体は
  他種と同様に動く (不活性マーカー) が、格子効果は一切モデルしない
- 酸化膜/窒化膜上の poly 核生成 (throw で拒否; poly が欲しければ
  deposit("poly", ...) を使う)
- オートドーピング (気相経由の基板→エピ輸送)・遷移層幅モデル —
  界面遷移は anneal の拡散のみで表現
- 連続する Si 層のタグ統合 — 各 epitaxy 呼び出しは独自タグを保持
  (テスト 3 が積層順の検証にタグ差を使う)
- oxidize 型のサブステップ分割 — アニールは末尾 1 回 (確定方針)
