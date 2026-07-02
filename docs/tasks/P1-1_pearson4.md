# P1-1: Pearson-IV 解析注入 (深さプロファイル)

## 目的

解析注入 (src/implant.cpp) は縦方向 Gaussian のみで、実測プロファイルの
非対称性 (表面側スキュー) と裾の重さを表現できない。モーメントテーブルを
4 モーメント {Rp, ΔRp, γ (歪度), β (尖度)} に拡張し、Pearson-IV 深さ
プロファイルを実装する。tilt/rotation と dual-Pearson (チャネリングテール)
は本タスクの対象外 — 深さプロファイルのみ。

## 現状コード

- `include/cprocess/materials.hpp` — `Dopant::range` は
  `std::vector<std::array<double, 3>>` ({E keV, Rp cm, dRp cm});
  `implant_range(d, E, rp, drp)` が log-E 線形補間 (範囲外はクランプ)
- `src/materials.cpp` — `kDopants` テーブル (B/P/As: 8 点、Sb: 5 点、
  NM=1e-7 で nm→cm)
- `include/cprocess/implant.hpp` — `ImplantParams { dopant, dose, rp, drp,
  drl, has_window, x1..y2 }`; `apply_implant(mesh, mask, p, conc)` が
  深さ d = ztop − z の Gaussian × 窓の erf 横方向係数を加算し原子数を返す
- `src/implant.cpp` — `apply_implant` 実装 (ピーク = dose/(√(2π)·dRp))
- `src/process.cpp` — `proc::implant_gauss(st, species, dose, energy_kev,
  rp, drp, drl, has_window, x1..y2, seed_damage, log)`:
  energy_kev>0 なら `implant_range` でテーブル引き、でなければ rp/drp 手動
- `src/deck.cpp` — implant コマンド (method=mc でなければ implant_gauss)
- `python/_cprocess.cpp` — `ImplantParams` バインディング、
  `proc_implant_gauss`
- `python/cprocess/simulation.py` — `Simulation.implant(..., mc=False, ...)`

## 実装手順

1. **モーメントテーブル拡張** (`include/cprocess/materials.hpp`):
   `Dopant::range` を `std::vector<std::array<double, 5>>` に変更
   ({E keV, Rp cm, dRp cm, gamma, beta})。`src/materials.cpp` の
   `kDopants` の各エントリに γ/β を追記する。値は教科書的な典型値
   (軽イオンほど・高エネルギーほど後方散乱で表面側スキューが強い。
   精度より単調傾向が重要):

   | E [keV] | B γ | B β | P γ | P β | As γ | As β |
   |---|---|---|---|---|---|---|
   | 10  | -0.5 | 3.5 | -0.30 | 3.2 | -0.20 | 3.1 |
   | 20  | -0.7 | 4.0 | -0.40 | 3.4 | -0.25 | 3.2 |
   | 30  | -0.8 | 4.5 | -0.45 | 3.5 | -0.30 | 3.3 |
   | 50  | -1.0 | 5.0 | -0.55 | 3.7 | -0.35 | 3.4 |
   | 80  | -1.2 | 6.0 | -0.65 | 4.0 | -0.45 | 3.6 |
   | 100 | -1.3 | 6.5 | -0.70 | 4.2 | -0.50 | 3.7 |
   | 150 | -1.4 | 7.5 | -0.80 | 4.6 | -0.60 | 4.0 |
   | 200 | -1.5 | 8.0 | -0.90 | 5.0 | -0.70 | 4.3 |

   Sb (5 点: 10/30/50/100/200 keV):
   γ = -0.20 / -0.30 / -0.35 / -0.50 / -0.60、
   β = 3.1 / 3.3 / 3.4 / 3.7 / 4.0。
   (全点で β > 1 + γ² を満たすことを確認済み: 最悪 B 200 keV で
   1+2.25=3.25 < 8.0)
2. **モーメント補間** (`materials.hpp` / `materials.cpp`):
   ```cpp
   // 4-moment lookup, log-E linear interpolation (clamped like implant_range).
   bool implant_moments(const Dopant& d, double energy_kev,
                        double& rp, double& drp, double& gamma, double& beta);
   ```
   実装は既存 `implant_range` と同じ補間ロジックで 4 成分を補間する。
   `implant_range` は `implant_moments` を呼んで rp/drp だけ返すラッパに
   書き換え (重複排除)。既存の `t[i][1]`/`t[i][2]` 参照はそのまま有効
3. **ImplantParams 拡張** (`include/cprocess/implant.hpp`):
   ```cpp
   enum class Profile { gauss, pearson4 };
   struct ImplantParams {
     ...(既存)...
     Profile profile = Profile::gauss;
     double gamma = 0.0, beta = 3.0;  // Pearson-IV moments (profile=pearson4)
   };
   ```
4. **Pearson-IV プロファイル** (`src/implant.cpp`, anonymous namespace)。
   Pearson 系のパラメータ (これをそのまま実装する; 導出不要):
   ```
   A  = 10*beta - 12*gamma^2 - 18
   b0 = -drp*drp * (4*beta - 3*gamma^2) / A
   b1 = -gamma * drp * (beta + 3) / A
   b2 = -(2*beta - 3*gamma^2 - 6) / A
   ```
   **有効条件** (満たさなければ Gaussian にフォールバック):
   `beta > 1 + gamma*gamma` かつ `b1*b1 - 4*b0*b2 < 0`
   (後者は分母 b0 + b1 z' + b2 z'^2 が実根を持たない Type-IV 条件)。
   プロファイルは f'(z')/f(z') = −(z' + b1)/(b0 + b1 z' + b2 z'^2)
   (z' = z − Rp) を数値積分した表として構築する:
   ```cpp
   struct Pearson4Table {
     double d_lo = 0, dz = 0;        // depth grid start / spacing [cm]
     std::vector<double> f;          // n points, integral(f * dz) = 1 (cm^-1)
     double eval(double depth) const;  // linear interp; 0 outside the grid
   };
   // Returns false when the moments are invalid (caller falls back to Gauss).
   bool build_pearson4(double rp, double drp, double gamma, double beta,
                       Pearson4Table& out);
   ```
   `build_pearson4` の手順:
   1. 上記 A, b0, b1, b2 を計算し、有効条件を検査 (不成立 → return false)
   2. 深さ格子: `d_lo = std::max(0.0, rp - 6*drp)`, `d_hi = rp + 6*drp`,
      n = 2000 点等間隔、`dz = (d_hi - d_lo)/(n-1)`
   3. ln f を台形則で構築: 格子点 i の z' = (d_lo + i*dz) − rp、
      `g(z') = -(z' + b1)/(b0 + b1*z' + b2*z'*z')` として、
      Rp に最も近い格子点を基準 (ln f = 0) に前後へ
      `lnf[i±1] = lnf[i] ± 0.5*dz*(g(z'_i) + g(z'_{i±1}))` と積分
   4. `f[i] = exp(lnf[i])`、台形則で `S = Σ` を計算し `f[i] /= (S*dz)` と
      正規化 (格子上の積分 = 1)。これにより表面 (d<0) 側の切り落とし分も
      込みでドーズが厳密に保存される
   5. `eval(depth)`: depth が [d_lo, d_hi] 外なら 0、内なら線形補間
5. **apply_implant の分岐** (`src/implant.cpp`):
   `p.profile == Profile::pearson4` のとき `build_pearson4(p.rp, p.drp,
   p.gamma, p.beta, tbl)` を試み、成功なら縦方向係数を
   `v = p.dose * tbl.eval(d)` に置き換える (失敗なら従来 Gaussian)。
   横方向の窓 erf 係数・mask 処理・atoms 積算は共通のまま変更しない
6. **proc:: 拡張** (`include/cprocess/process.hpp` / `src/process.cpp`):
   `implant_gauss` に引数を 1 つ追加 (seed_damage の後):
   ```cpp
   double implant_gauss(SimState& st, const std::string& species, double dose,
                        double energy_kev, double rp, double drp, double drl,
                        bool has_window, double x1, double x2, double y1,
                        double y2, bool seed_damage = false,
                        const std::string& profile = "gauss",
                        std::ostream* log = nullptr);
   ```
   - `profile` は "gauss" | "pearson" (それ以外は
     `std::runtime_error("implant: profile must be gauss|pearson")`)
   - "pearson" は `energy_kev > 0` を必須とし (rp/drp 手動指定と併用不可 →
     throw "pearson profile requires energy=")、`implant_moments` で
     rp/drp/gamma/beta を取得して `p.profile = Profile::pearson4` を設定
   - フォールバック発生の可視化のため、ログ行に `profile=pearson4` または
     `profile=gauss` を出す
   - `src/deck.cpp` の implant コマンドに `profile=` キー (既定 "gauss") を
     追加し、そのまま proc へ渡す
7. **pybind11** (`python/_cprocess.cpp`):
   - `proc_implant_gauss` のラムダとシグネチャに
     `py::arg("profile") = "gauss"` を追加 (damage の後)
   - `ImplantParams` バインディングに追加:
     ```cpp
     .def_readwrite("gamma", &ImplantParams::gamma)
     .def_readwrite("beta",  &ImplantParams::beta)
     .def_property("profile",
         [](const ImplantParams& p) {
           return p.profile == ImplantParams::Profile::pearson4
                      ? "pearson" : "gauss"; },
         [](ImplantParams& p, const std::string& s) {
           if (s == "pearson") p.profile = ImplantParams::Profile::pearson4;
           else if (s == "gauss") p.profile = ImplantParams::Profile::gauss;
           else throw std::runtime_error("profile must be gauss|pearson");
         })
     ```
8. **Simulation メソッド** (`python/cprocess/simulation.py`):
   `Simulation.implant` にキーワード `profile: str = "gauss"` を追加。
   - `mc=True and profile != "gauss"` → `ValueError("profile applies to the
     analytic implant only")`
   - 解析経路では `_c.proc_implant_gauss(..., bool(damage), profile)` と
     最後に渡す。docstring に「profile="pearson" は energy 指定必須、
     4 モーメント (Rp, ΔRp, γ, β) テーブルによる Pearson-IV 縦方向
     プロファイル」を追記
9. `tests/test_pearson.cpp` を新規作成し、`CMakeLists.txt` の
   `foreach(t ...)` に `pearson` を追加

## テスト仕様

`tests/test_pearson.cpp` (test_util.hpp 使用。細分メッシュ:
`proc::mesh_box(st, 0, 0.2e-4, 0, 0.2e-4, 0, 1.2e-4, 2, 2, 240)`;
深さモーメントはセル場から質量加重で計算:
mean = Σ cV·d / Σ cV, sigma² = Σ cV·d²/Σ cV − mean², d = ztop − z):

1. **ドーズ保存**: B 200 keV (Rp=531 nm, ΔRp=94 nm, γ=−1.5, β=8.0;
   Rp−6ΔRp = −33 nm なので表面切り落としは僅少), dose=1e14,
   profile="pearson" で `proc::implant_gauss`。戻り値 atoms と
   `dose * area` (area = 0.2e-4 × 0.2e-4 cm²) の相対差 **< 0.5%**。
   場の積分 Σ c_i V_i でも同じく < 0.5%
2. **モーメント一致**: 同じ注入で
   `|mean − 531e-7| < 0.02 * 531e-7` (**Rp ±2%**)、
   `|sigma − 94e-7| < 0.05 * 94e-7` (**ΔRp ±5%**)
3. **負の歪度**: 同プロファイルの 3 次中心モーメント
   μ3 = Σ cV·(d−mean)³ / Σ cV が **μ3 < 0** (裾が表面側) を CHECK。
   比較用に profile="gauss" の同注入では |skew| = |μ3|/σ³ < 0.1 を CHECK
4. **Gaussian フォールバック**: `apply_implant` を直接使い、
   `p.profile = Profile::pearson4, p.gamma = -2.0, p.beta = 4.0`
   (β=4.0 ≤ 1+γ²=5.0 → 無効) で例外なく走り、結果が
   `p.profile = Profile::gauss` の同条件実行と全セル bit 一致 (==) すること
5. **profile="pearson" + energy 無し (rp/drp 手動) は throw**、
   `profile="foo"` も throw
6. **既存 Gaussian の非退行**: profile 既定値での `proc::implant_gauss` が
   従来と bit 一致 (既存 test_diffusion / test_integration がそのまま PASS
   することで担保。テスト内でも profile 引数省略呼び出しを 1 回実行)

`python/test_comprehensive.py` に追加:

- `sim.implant("B", dose=1e14, energy=200, profile="pearson")` 後、
  `abs(sim.dose("B")/1e14 - 1) < 0.01` を assert
- 深さプロファイルの歪度が負 (`sim.field("B")`, `sim.cell_volumes`,
  `sim.cell_centroids[:,2]` から μ3 < 0) を assert
- `sim.implant("B", dose=1e13, rp=0.1, drp=0.03, profile="pearson")` が
  例外 (RuntimeError) になること、
  `sim.implant("B", dose=1e13, energy=50, mc=True, profile="pearson")` が
  `ValueError` になることを assert

実行: `cmake --build build -j$(nproc) && ctest --test-dir build` 全 PASS、
`build/_cprocess*.so` を `python/cprocess/` へコピーして
`python3 python/test_comprehensive.py`。

## 完了条件 (DoD)

- [ ] `Dopant::range` が array<double,5> になり、B/P/As/Sb 全点に上表の
      γ/β が入っている (`implant_moments` + `implant_range` ラッパ化)
- [ ] Pearson-IV: 上記 A/b0/b1/b2 閉形式 + Type-IV 有効条件 + 2000 点
      深さ表 + 正規化 (フォールバックは Gaussian、throw しない)
- [ ] `proc::implant_gauss` の `profile` 引数 + deck の `profile=` キー +
      pybind (`proc_implant_gauss` / `ImplantParams.profile/gamma/beta`) +
      `Simulation.implant(profile=...)`
- [ ] `tests/test_pearson.cpp` (上記 6 ケース) を CMakeLists の foreach に
      登録、既存テスト含め全 PASS
- [ ] Python テスト追加、`python3 python/test_comprehensive.py` PASS
- [ ] コミットメッセージに `P1-1` を含める

## やらないこと

- 解析注入の tilt/rotation 対応 — 深さプロファイルのみが本タスクの範囲
- dual-Pearson (チャネリングテール第 2 成分) — 将来タスク。チャネリングが
  必要なら MC (`mc=True, channeling=True`) を使う
- 横方向広がりモデルの変更 — 窓の erf 畳み込みは Gaussian のまま流用
- モーメントテーブルの精密較正・エネルギー点の追加 (P1-10 の上書き機構で
  対応予定)。値は単調傾向が正しければよい
- MC 注入・レジストスタック経路への影響 (ImplantParams は解析注入専用)
