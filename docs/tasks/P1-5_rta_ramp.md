# P1-5: RTA 温度ランプ (時間依存温度プロファイル)

## 目的

`diffuse` / `diffuse_ted` は等温のみ (`DiffuseOpts::temp` 一定) で、
RTA のスパイクアニール (例 900 °C → 1050 °C → 900 °C) が表現できない。
`DiffuseOpts` に区分線形の温度プロファイル T(t) を追加し、各時間ステップで
中点温度 T(t_mid) を評価して ni・拡散係数・(run_ted では) 格子間原子
パラメータを再計算する。時間刻みはプロファイルの折れ点をまたがないよう
クランプする。

## 現状コード

- `include/cprocess/diffusion.hpp` — `DiffuseOpts { temp, time, dt, ... }`。
  dt<=0 は time/50
- `src/diffusion.cpp` — `DiffusionSolver::run()`:
  - ループ前に一度だけ `const double dt0 = ...; const int nsteps =
    static_cast<int>(std::ceil(o.time / dt0 - 1e-12)); const double ni =
    ni_si(o.temp);` を計算
  - 時間ループは `for (int step = 0; step < nsteps; ++step) { const double
    dt = std::min(dt0, o.time - t); ... }` の固定ステップ数
  - Picard 内で `dopant_diffusivity(dp, o.temp, nni[i])` が温度を直接参照
- `src/diffusion.cpp` — `DiffusionSolver::run_ted()`: 同様の構造に加え、
  ループ前に `cstar = interstitial_cstar(o.temp)`, `d_I =
  interstitial_diffusivity(o.temp)`, `k_rec = interstitial_recomb_rate(o.temp)`
  と格子間原子拡散係数ベクトル `dI` (mask 依存) を一度だけ計算している。
  ループ内では `assemble(dI, ..., k_rec, ...)`、S の計算に `cstar` を使用
- `src/process.cpp` — `proc::diffuse` / `proc::diffuse_ted` は opts を
  そのまま渡し、最後に `st.last_temp = opts.temp;`
- `python/_cprocess.cpp` — `DiffuseOpts` バインディング (temp/time/dt/...)
- `python/cprocess/simulation.py` — `Simulation.diffuse(time, temp, *, dt,
  field_enh, nonortho, ted)`; time/dt は分、temp は °C

## 実装手順

1. `include/cprocess/diffusion.hpp` の `DiffuseOpts` に追加:
   ```cpp
   // Piecewise-linear temperature profile {time_s from step start, temp_K}.
   // Empty = isothermal at `temp`. Must be sorted, start at t=0, size >= 2.
   // Beyond the last breakpoint the last temperature is held.
   std::vector<std::pair<double, double>> temp_profile;
   ```
   同ヘッダに自由関数を宣言 (process.cpp からも使うため public):
   ```cpp
   // T(t) [K]: linear interpolation of opts.temp_profile; opts.temp if empty.
   double temp_at(const DiffuseOpts& opts, double t);
   ```
2. `src/diffusion.cpp` に `temp_at` を実装:
   - `o.temp_profile.empty()` → `return o.temp;`
   - t <= 最初の点 → 最初の T、t >= 最後の点 → 最後の T、それ以外は
     隣接 2 点の線形補間
   さらにプロファイル検証ヘルパー (anonymous namespace)
   `void validate_profile(const DiffuseOpts& o)`:
   非空なら size>=2・`temp_profile[0].first == 0`・時刻が狭義単調増加・
   全温度 > 0 を確認、違反時 `std::runtime_error("diffuse: bad temp_profile")`。
   `run()` / `run_ted()` の冒頭で呼ぶ。
3. `run()` の時間ループを可変ステップ化する:
   - `const int nsteps = ...` の行を削除。`const double dt0` は現状どおり
   - 折れ点クランプ付き while ループに置換:
     ```cpp
     double t = 0;
     int step = 0;
     while (t < o.time - 1e-12 * o.time) {
       double dt = std::min(dt0, o.time - t);
       for (const auto& [tb, Tb] : o.temp_profile)   // next breakpoint clamp
         if (tb > t + 1e-12 * o.time && tb - t < dt) dt = tb - t;
       const double T = temp_at(o, t + 0.5 * dt);    // mid-step temperature
       const double ni = ni_si(T);
       ...   // (existing per-step body)
       t += dt;
       ++step;
     }
     ```
   - ループ前の `const double ni = ni_si(o.temp);` を削除し、上記のとおり
     ループ内へ移動。Picard 内の `dopant_diffusivity(dp, o.temp, nni[i])` を
     `dopant_diffusivity(dp, T, nni[i])` に変更 (field_enh ブロックは nni のみ
     依存なので変更不要。P1-3 実装済みなら activation クランプも
     `active_concentration(..., T)` に変更)
   - ログ: 総ステップ数が事前に確定しないため、`step %4d/%d` を
     `step %4d  t=%.6g s  T=%.5g K  picard=%d  cg=%d` 形式に変更。
     出力頻度は `verbosity >= 2` なら毎ステップ、`verbosity == 1` なら
     従来同様おおよそ 10 行 (`every = std::max(1, (int)std::ceil(o.time/dt0) / 10)`)
4. `run_ted()` も同じ while ループ構造に変更し、さらに次の 4 つを
   ループ内 (ステップ冒頭、dt と T の確定直後) で再計算する:
   ```cpp
   const double ni    = ni_si(T);
   const double cstar = interstitial_cstar(T);
   const double d_I   = interstitial_diffusivity(T);
   const double k_rec = interstitial_recomb_rate(T);
   for (int i = 0; i < nc; ++i) dI[i] = mask_[i] ? d_I : 0.0;
   ```
   (`dI` ベクトル自体はループ外で確保したまま、値の更新のみループ内)。
   `dopant_diffusivity(dp, o.temp, ...)` → `(dp, T, ...)`。
   ループ後のサマリログの `cstar` 参照は最後のステップの値でよい
5. `proc::diffuse` / `proc::diffuse_ted` (src/process.cpp):
   `st.last_temp = opts.temp;` を
   `st.last_temp = temp_at(opts, opts.time);` に変更 (終了時温度)。
   先頭ログにプロファイル使用時は
   `T=ramp(900->1050->900 K)` のように折れ点温度列を出す (任意の簡易書式で可)
6. pybind11 (`python/_cprocess.cpp`): `DiffuseOpts` バインディングに
   ```cpp
   .def_readwrite("temp_profile", &DiffuseOpts::temp_profile,
       "Piecewise-linear T profile [(t_s, T_K), ...]; empty = isothermal")
   ```
   を追加 (`pybind11/stl.h` は既に include 済みで
   `std::vector<std::pair<double,double>>` は自動変換される)。
   proc 関数のシグネチャは変更なし (opts が運ぶ)
7. `python/cprocess/simulation.py` — `Simulation.diffuse` を拡張:
   ```python
   def diffuse(self, time: float, temp: float = None, *,
               dt: float = 0.0, field_enh: bool = True,
               nonortho: bool = True, ted: bool = False,
               ramp=None) -> "Simulation":
   ```
   - `ramp` は `[(t_min, T_celsius), ...]`。指定時は
     `ramp[0][0] == 0` を要求 (違反は `ValueError("ramp must start at t=0")`)、
     `opts.temp_profile = [(t * MIN, _celsius_to_k(T)) for t, T in ramp]`、
     `opts.temp = _celsius_to_k(ramp[0][1])` を設定し、**引数 temp は無視**する
     (docstring に明記)
   - `ramp is None and temp is None` → `ValueError("give temp or ramp")`
   - ramp なしの挙動は従来どおり
8. `tests/test_rta.cpp` を新規作成し、`CMakeLists.txt` の `foreach(t ...)` に
   `rta` を追加

## テスト仕様

`tests/test_rta.cpp` (tests/test_ted.cpp の `profile_spread` ヘルパーを
コピーして使用。共通セットアップ: `proc::mesh_box(st, 0, 0.3e-4, 0, 0.3e-4,
0, 1.0e-4, 4, 4, 40)` + B Gaussian 注入 `dose=1e14, rp=0.05e-4, drp=0.02e-4`):

1. **定数プロファイル = 等温と一致**: `opts.time=60, opts.temp=1273.15,
   verbosity=0` の等温 run と、同一初期状態で
   `opts.temp_profile={{0,1273.15},{60,1273.15}}` の run を比較。
   折れ点 (0 s, 60 s) はステップ境界と一致するので刻みは同一になり、
   全セルで**相対差 < 1e-12** (`|a-b| <= 1e-12*(|a|+1e-300)`) を CHECK
2. **スパイクは中間**: `time=120` s で 3 通り:
   (a) 等温 900 °C (`temp=1173.15`)、(b) 等温 1050 °C (`temp=1323.15`)、
   (c) ランプ `{{0,1173.15},{60,1323.15},{120,1173.15}}`。
   拡散後 spread が **spread(a) < spread(c) < spread(b)** (狭義) を CHECK
3. **折れ点クランプ**: `time=10, dt=3.0`、
   `temp_profile={{0,1173.15},{1.0,1273.15},{10,1273.15}}`、`verbosity=2` で
   ログを `std::ostringstream` に取り、`[diffuse]   step` を含む行数が
   **ちょうど 4** (ステップ終了時刻 1, 4, 7, 10 s) であること、および
   1 行目に `t=1 ` (t=1.0 s ちょうど) が現れることを CHECK
4. **run_ted のランプ**: damage=true の B 注入
   (`proc::implant_gauss(..., seed_damage=true)`) 後、
   `proc::diffuse_ted` をランプ `{{0,1173.15},{30,1323.15},{60,1173.15}}`,
   `time=60` で実行。例外なく完走し、質量保存 (B 総原子数の相対変化 < 1e-6)、
   かつ spread が等温 900 °C 60 s の diffuse_ted より大きいことを CHECK
5. **不正プロファイルは throw**: `temp_profile={{5,1273.15},{10,1273.15}}`
   (t=0 で始まらない) と `{{0,1273.15}}` (1 点) がそれぞれ
   `std::runtime_error` になることを CHECK
6. **last_temp**: テスト 2(c) 実行後 `st.last_temp == 1173.15` (終了時温度)

`python/test_comprehensive.py` に追加:

- `sim.diffuse(time=2, ramp=[(0, 900), (1, 1050), (2, 900)])` が完走し
  dose が保存される (相対変化 < 1e-4) こと
- `sim.diffuse(time=2, ramp=[(0.5, 900), (2, 900)])` が `ValueError`、
  `sim.diffuse(time=2)` (temp も ramp もなし) が `ValueError` になること
- `sim.diffuse(time=1, temp=1000)` (従来 API) が引き続き動くこと

実行: `cmake --build build -j$(nproc) && ctest --test-dir build` 全 PASS、
`build/_cprocess*.so` を `python/cprocess/` へコピーして
`python3 python/test_comprehensive.py`。

## 完了条件 (DoD)

- [ ] `DiffuseOpts::temp_profile` + `temp_at()` + 検証 (throw)
- [ ] `run()` / `run_ted()` の**両方**が可変ステップ while ループになり、
      ni / 拡散係数 / (run_ted の) cstar・d_I・k_rec・dI を毎ステップ再計算
- [ ] 折れ点をまたがない dt クランプ (テスト 3 で検証)
- [ ] `proc::diffuse` / `diffuse_ted` の `last_temp` が終了時温度になる
- [ ] pybind `temp_profile` + `Simulation.diffuse(ramp=...)` (ramp[0][0]==0
      必須、temp 無視) + Python テスト
- [ ] `tests/test_rta.cpp` (上記 6 ケース) を CMakeLists の foreach に登録、
      全テスト PASS
- [ ] コミットメッセージに `P1-5` を含める

## やらないこと

- 適応時間刻み (局所誤差推定) — S-5。本タスクは折れ点クランプのみ
- ランプレート指定 (°C/s) の糖衣 API — Python 側は折れ点リストのみ
- デッキ (deck.cpp) への ramp 構文追加 — Python 前面で足りる。必要なら別タスク
- 温度依存の固溶度クランプとの結合調整 — P1-3 が o.temp→T を差し替えるだけ
  (両タスクどちらが先でも可; 後から入れる側が 1 行合わせる)
- 空間的に不均一な温度場
