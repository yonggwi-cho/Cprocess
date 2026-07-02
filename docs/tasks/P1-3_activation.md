# P1-3: 活性化モデル (拡散中の固溶度クランプ)

## 目的

現状、固溶度クランプは `proc::save` (src/process.cpp) が VTK 出力用に
`X_active = min(C, C_ss(last_temp))` を計算するときにしか使われず、
**拡散計算そのものは全濃度 C を電気的に活性として扱う**。高ドーズ As/P では
ピーク濃度が固溶度を大きく超えるため、電荷中性から求める n/ni と電場増速が
過大になり、拡散が速すぎる。SUPREM 流の最も単純な活性化モデル
「拡散中も活性濃度 C_act = min(C, C_ss(T)) だけが電荷中性に寄与する」を導入する。
全濃度 C は引き続き拡散する (不活性分を不動にするクラスタ動力学は P2-2)。

## 現状コード

- `src/materials.cpp` — `solid_solubility(const Dopant&, double temp_k)`:
  `ss_pre * exp(-ss_e / kT)`。ss_pre<=0 なら 0 を返す (0 = 「固溶度なし」)
- `src/diffusion.cpp` — `DiffusionSolver::run()` の Picard ループ冒頭に
  電荷中性ループがある:
  ```cpp
  for (int i = 0; i < nc; ++i) {
    double nnet = 0;
    for (int s = 0; s < ns; ++s) {
      const double c = (*fields[s].conc)[i];
      nnet += (fields[s].dopant->type == DopType::donor) ? c : -c;
    }
    const double cc = nnet / (2.0 * ni);
    nni[i] = cc + std::sqrt(cc * cc + 1.0);
  }
  ```
  **同一のループが `run_ted()` にもある** (Picard ループ内、`for (int picard = 1; ...)`
  直後)。両方を変更する。nni は続けて `dopant_diffusivity(dp, o.temp, nni[i])` と
  電場増速係数 (`o.field_enh` ブロック) に使われる — これらは nni 経由で自動的に
  クランプの影響を受けるので追加変更は不要
- `src/process.cpp` — `proc::save()` が各種の `X_active` を
  `std::min(conc[i], css)` で計算 (css = `solid_solubility(*d, st.last_temp)`)。
  この計算を新ヘルパーに置き換えて一元化する。deck.cpp / vtk_writer.cpp に
  他のクランプ箇所はない
- `include/cprocess/diffusion.hpp` — `DiffuseOpts` (temp/time/dt/field_enh/...)
- `SimState::last_temp` (include/cprocess/deck.hpp) — 最後の diffuse 温度 [K]、
  初期値 1273.15

参考数値 — 既存の ss_pre/ss_e 定数 (src/materials.cpp) から計算した固溶度:

| 種 | ss_pre [cm^-3] | ss_e [eV] | C_ss(900 °C) | C_ss(1000 °C) |
|----|---------------|-----------|--------------|---------------|
| B  | 9.25e22 | 0.73 | 6.7e19 | 1.19e20 |
| P  | 2.45e23 | 0.62 | 5.3e20 | 8.6e20 |
| As | 1.3e23  | 0.66 | 1.90e20 | 3.17e20 |
| Sb | 3.8e21  | 0.56 | 1.49e19 | 2.31e19 |

(kT(900 °C=1173.15 K)=0.10110 eV, kT(1000 °C=1273.15 K)=0.10971 eV)

## 実装手順

1. `include/cprocess/materials.hpp` に宣言、`src/materials.cpp` に実装:
   ```cpp
   // Electrically active concentration: solid-solubility clamp
   //   C_act = min(C, C_ss(T));  C_ss == 0 (no fit) means "no clamp".
   double active_concentration(const Dopant& d, double conc, double temp_k);
   ```
   実装は `const double css = solid_solubility(d, temp_k);
   return (css > 0) ? std::min(conc, css) : conc;`
2. `include/cprocess/diffusion.hpp` の `DiffuseOpts` にフラグを追加:
   ```cpp
   bool activation = true;  // clamp charge neutrality at solid solubility
   ```
   (テスト・回帰比較用に off にできるようにする)
3. `src/diffusion.cpp` — `run()` と `run_ted()` の**両方**の電荷中性ループで、
   `const double c = (*fields[s].conc)[i];` を次に置き換える:
   ```cpp
   double c = (*fields[s].conc)[i];
   if (o.activation) c = active_concentration(*fields[s].dopant, c, o.temp);
   ```
   それ以外 (nnet 加算、nni 式、dcell ループ、field_enh) は変更しない。
   ※ P1-5 (RTA ランプ) 実装後は `o.temp` が「そのステップの T(t)」に
   変わるが、本タスクでは `o.temp` のままでよい (両タスクは独立)。
4. `src/process.cpp` — `proc::save()` の
   `act[i] = (css > 0) ? std::min(conc[i], css) : conc[i];` を
   `act[i] = active_concentration(*d, conc[i], st.last_temp);` に置き換え
   (css のローカル計算は削除)。
5. 活性濃度クエリを `proc::` に追加。`include/cprocess/process.hpp` に宣言、
   `src/process.cpp` に実装:
   ```cpp
   // Per-cell electrically active concentration of `species` at temp_k [K].
   // temp_k <= 0 selects st.last_temp. Throws if the field does not exist.
   std::vector<double> active_field(const SimState& st,
                                    const std::string& species,
                                    double temp_k = -1.0,
                                    std::ostream* log = nullptr);
   ```
   実装: `dopant_or_throw(species)` → `st.fields.at(d->symbol)` (無ければ
   `std::runtime_error("no field: " + species)`) → 各セルに
   `active_concentration` を適用した vector を返す。log には
   `[active] X: C_ss(T)=... cm^-3` を 1 行出す。
6. pybind11 バインディング (`python/_cprocess.cpp`):
   - `DiffuseOpts` バインディングに
     `.def_readwrite("activation", &DiffuseOpts::activation)` を追加
   - クエリ関数 (`proc_*` 規約。クエリなので log ではなく配列を返す):
     ```cpp
     m.def("proc_active_field",
         [](const SimState& st, const std::string& species, double temp_k) {
           return vec_to_np(proc::active_field(st, species, temp_k, nullptr));
         },
         py::arg("state"), py::arg("species"), py::arg("temp_k") = -1.0);
     ```
7. `python/cprocess/simulation.py`:
   - `Simulation.diffuse(...)` にキーワード `activation: bool = True` を追加し
     `opts.activation = activation` を設定
   - クエリメソッドを追加 (データを返すので self ではなく配列を返す):
     ```python
     def active(self, species: str, temp: float = None) -> np.ndarray:
         """Electrically active concentration [cm^-3] (solid-solubility clamp).

         temp in Celsius; None uses the temperature of the last diffuse step.
         """
         t_k = -1.0 if temp is None else _celsius_to_k(temp)
         return _c.proc_active_field(self._st, species, t_k)
     ```
8. `tests/test_activation.cpp` を新規作成し、`CMakeLists.txt` の
   `foreach(t ...)` リストに `activation` を追加する。

## テスト仕様

`tests/test_activation.cpp` (test_util.hpp の CHECK を使用):

1. **クランプで拡散が遅くなる**: 箱メッシュ
   `proc::mesh_box(st, 0, 0.3e-4, 0, 0.3e-4, 0, 1.0e-4, 4, 4, 40)`、
   高ドーズ As Gaussian 注入
   `proc::implant_gauss(st, "As", 3e15, 0, 0.05e-4, 0.02e-4, ...)`
   (ピーク = 3e15/(√(2π)·2e-6) ≈ 6.0e20 cm^-3 > C_ss(900 °C)=1.90e20)。
   同一初期状態 2 つを `DiffuseOpts{temp=1173.15, time=1800 (30 min),
   verbosity=0}` で、片方 `activation=true`、片方 `false` で
   `proc::diffuse`。tests/test_ted.cpp の `profile_spread` と同じ
   質量加重 z 標準偏差で比較し、
   **`spread_off > 1.02 * spread_on`** (2% 超の差) を CHECK。
   両者とも質量保存 (総原子数の相対変化 < 1e-6) を CHECK
2. **active ≤ field、低濃度側では等しい**: 上のクランプ有り状態で
   `auto act = proc::active_field(st, "As", 1173.15);` を取り、全セルで
   `act[i] <= C[i] + 1e-30` かつ `act[i] <= 1.90001e20`。さらに
   `C[i] < 1.8e20` のセルでは `act[i] == C[i]` (bit 一致) を CHECK
3. **低濃度なら active == field 全域一致**: 別状態で
   `proc::init(st, "B", 1e15)` のみ → `active_field(st, "B", 1273.15)` が
   全セルで field と bit 一致 (1e15 ≪ C_ss(1000 °C)=1.19e20)
4. **temp_k デフォルト**: `st.last_temp` を diffuse で設定した後、
   `active_field(st, "As")` (temp_k=-1) と
   `active_field(st, "As", st.last_temp)` が全セル bit 一致
5. **存在しない種は throw**: `active_field(st, "P")` (P 場なし) が
   `std::runtime_error` を投げることを try/catch で CHECK

`python/test_comprehensive.py` に追加:

- `sim.mesh(...)`; `sim.implant("As", dose=3e15, rp=0.05, drp=0.02)`;
  `sim.diffuse(time=30, temp=900)` 後に `a = sim.active("As")` として
  `np.all(a <= sim.field("As") + 1e-30)` かつ
  `a.max() <= 1.91e20` (C_ss(900 °C) 近傍で頭打ち) を assert
- `sim.active("As", temp=1000)` の最大値が `3.18e20` 以下かつ
  `1.9e20` 超 (温度依存が効いている) を assert
- `sim.diffuse(..., activation=False)` が例外なく走ることを確認

実行: `cmake --build build -j$(nproc) && ctest --test-dir build` 全 PASS、
`build/_cprocess*.so` を `python/cprocess/` にコピーして
`python3 python/test_comprehensive.py`。

## 完了条件 (DoD)

- [ ] `active_concentration` が materials に追加され、`proc::save` も
      それを使う (クランプ実装の一元化)
- [ ] `run()` / `run_ted()` の**両方**の電荷中性ループがクランプ値を使う
      (`DiffuseOpts::activation` で切替可、既定 on)
- [ ] `proc::active_field` + `proc_active_field` バインディング +
      `Simulation.active()` (クエリなので配列を返す)
- [ ] `tests/test_activation.cpp` (上記 5 ケース) を CMakeLists の foreach に
      登録、全テスト PASS
- [ ] Python テスト追加、`python3 python/test_comprehensive.py` PASS
- [ ] コミットメッセージに `P1-3` を含める

## やらないこと

- クラスタ動力学 (BIC, As-V) と活性化の時間依存 — P2-2。本タスクの
  クランプは瞬時平衡
- 不活性分を不動にする (拡散フラックスから除く) こと — 全濃度 C が
  従来どおり拡散する。dcell の変更は nni 経由の間接効果のみ
- `dopant_diffusivity` 自体や ss_pre/ss_e 定数の再較正
- デッキ (deck.cpp) への新コマンド追加 (activation は既定 on で十分)
- 温度ランプとの結合 (P1-5 が独立に実装する)
