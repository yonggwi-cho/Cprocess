# P1-10: パラメータ上書き機構 (ParamDB, pdbSet 相当)

## 目的

物性値 (拡散係数プレファクタ、固溶度、点欠陥パラメータ等) がすべて
`src/materials.cpp` にハードコードされており、較正・感度解析ができない。
実行時に上書き可能な**パラメータレジストリ `cp::ParamDB`** を導入し、
materials.cpp の各関数が「上書き値があればそれを、なければコンパイル時定数を」
読むようにする。Python から `sim.set_param("B.d0", ...)` で較正できることが
ゴール (SProcess の pdbSet 相当の最小版)。

## 現状コード

- `src/materials.cpp` — 無名 namespace の static テーブル `kDopants`
  (B/P/As/Sb の d0/e0, dm/em, dmm/emm, dp/ep, ss_pre/ss_e, fi, range)。
  `dopant_diffusivity(d, T, n/ni)`, `solid_solubility(d, T)`,
  `interstitial_cstar(T)` (3.0e27, 3.7 eV),
  `interstitial_diffusivity(T)` (5.0e-2, 1.77 eV),
  `interstitial_recomb_rate(T)` (2.0e4, 1.4 eV) は定数を直接埋め込み
- `include/cprocess/materials.hpp` — `Dopant` struct (symbol フィールドあり)、
  上記関数宣言
- TED 定数: `smax` 上限・Frenkel 生存率相当の定数は `src/diffusion.cpp`
  の run_ted 内にある — grep で `smax` / survival 相当の定数を特定して
  同様に ParamDB 経由に置き換える (キー名は下表)
- `python/_cprocess.cpp` — `PYBIND11_MODULE(_cprocess, m)` フラット構成、
  `proc_*` 関数群
- `python/cprocess/simulation.py` — `Simulation` クラス

## 実装手順

1. 新規 `include/cprocess/param_db.hpp` + `src/param_db.cpp`:
   ```cpp
   namespace cp {
   // 実行時パラメータ上書きレジストリ (シングルトン)。
   // スレッド安全性: set は solve 実行中に呼んではならない (ロックなし。
   // 制御フローは単一スレッドであり、solve 内の OpenMP 並列区間は
   // get を read-only にしか使わない)。この制約をヘッダに文書化する。
   class ParamDB {
    public:
     static ParamDB& instance();
     // 上書きがあればその値、なければ fallback を返す。
     double get(const std::string& key, double fallback) const;
     void set(const std::string& key, double value);
     bool erase(const std::string& key);       // 上書き解除
     void clear();                              // 全解除 (テスト用)
     std::map<std::string, double> all() const; // 現在の上書き一覧
    private:
     ParamDB() = default;
     std::map<std::string, double> overrides_;
   };
   }  // namespace cp
   ```
   CMakeLists.txt の `cprocess_core` ソースに `src/param_db.cpp` を追加。
2. **キー命名規約** (これを正とし、materials.hpp のコメントにも記載):

   | キー | 消費箇所 | fallback |
   |------|---------|----------|
   | `<Sym>.d0` / `<Sym>.e0` | dopant_diffusivity の中性項 | d.d0 / d.e0 |
   | `<Sym>.dm` / `<Sym>.em` | 同・単負項 | d.dm / d.em |
   | `<Sym>.dmm` / `<Sym>.emm` | 同・二重負項 | d.dmm / d.emm |
   | `<Sym>.dp` / `<Sym>.ep` | 同・正項 | d.dp / d.ep |
   | `<Sym>.ss_pre` / `<Sym>.ss_e` | solid_solubility | d.ss_pre / d.ss_e |
   | `<Sym>.fi` | run_ted の (1−fi)+fi·S | d.fi |
   | `I.cstar_pre` / `I.cstar_e` | interstitial_cstar | 3.0e27 / 3.7 |
   | `I.d0` / `I.e0` | interstitial_diffusivity | 5.0e-2 / 1.77 |
   | `I.krec_pre` / `I.krec_e` | interstitial_recomb_rate | 2.0e4 / 1.4 |
   | `ted.smax` | run_ted の過飽和上限 S_max | 現行定数 |
   | `ted.frenkel_survival` | "+1" シード時の生存率係数 | 現行値 (係数が明示されていなければ 1.0 を fallback とし seed_interstitials に乗算を追加) |

   `<Sym>` は `Dopant::symbol` そのまま ("B", "P", "As", "Sb")。
   単位はすべて**コア単位** (cm²/s, eV, cm⁻³, 1/s, 無次元) — 変換しない。
3. **materials.cpp の改修** (kDopants は再構築**しない**。上書きは
   読み出し時にレイヤする):
   - `dopant_diffusivity`: 各項の使用直前に
     `const auto& P = ParamDB::instance();
      const double d0 = P.get(d.symbol + ".d0", d.d0);` の形で 10 値
     (d0/e0/dm/em/dmm/emm/dp/ep) を取得して従来式に代入。
     **注意**: 現行は `if (d.d0 > 0)` でガードしている — ガードも
     上書き後の値で判定する (`if (d0 > 0)`) こと。0 → 正の上書きで項を
     有効化できるようにする
   - `solid_solubility`: ss_pre/ss_e を同様に。`ss_pre <= 0` ガードも
     上書き後の値で判定
   - `interstitial_cstar` / `interstitial_diffusivity` /
     `interstitial_recomb_rate`: prefactor と活性化エネルギーを
     `P.get("I.cstar_pre", 3.0e27)` 等に置換
   - `find_dopant` / `dopant_table` は**無変更** (const 参照を返したまま。
     copy-on-read はしない — fi の消費箇所 run_ted 側で
     `P.get(d->symbol + ".fi", d->fi)` を読む)
   - `src/diffusion.cpp` run_ted: fi と smax の読み出しを ParamDB 経由に
     置換 (時間ループの**外**で 1 回読むこと — ループ内 get は避ける)。
     `src/process.cpp` seed_interstitials: 加算量に
     `P.get("ted.frenkel_survival", 1.0)` を乗じる
4. **proc:: API** (process.hpp / process.cpp):
   ```cpp
   void set_param(SimState& st, const std::string& key, double value,
                  std::ostream* log = nullptr);   // st は未使用 (規約統一のため受ける)
   double get_param(const std::string& key, double fallback);
   std::map<std::string, double> list_params();   // 上書き中の一覧
   ```
   set_param はログに `[param] <key> = <value>` を出力。キーの妥当性検査は
   **しない** (未知キーも保持 — 消費側が読まなければ無効なだけ。
   docstring に明記)。
5. **pybind** (python/_cprocess.cpp):
   ```cpp
   m.def("proc_set_param", [](SimState& st, const std::string& k, double v) {
       std::ostringstream log; proc::set_param(st, k, v, &log); return log.str(); },
     py::arg("state"), py::arg("key"), py::arg("value"));
   m.def("proc_get_param", [](const std::string& k, double fb) {
       return proc::get_param(k, fb); },
     py::arg("key"), py::arg("fallback") = 0.0);
   m.def("proc_list_params", []() { return proc::list_params(); });
   ```
6. **Simulation** (python/cprocess/simulation.py):
   ```python
   def set_param(self, key: str, value: float) -> "Simulation":
       """Override a physical parameter (raw core units: cm^2/s, eV, cm^-3).
       No unit conversion is applied. e.g. set_param("B.d0", 0.074)."""
       self._emit(_c.proc_set_param(self._st, key, float(value)))
       return self

   def get_param(self, key: str, fallback: float = 0.0) -> float:
       return _c.proc_get_param(key, float(fallback))
   ```
   (set_param はチェーン、get_param はクエリなので値を返す — CLAUDE.md 規約)
7. テスト: 新規 `tests/test_params.cpp` を作成し CMakeLists.txt の
   foreach に `params` を追加。各テストの末尾で
   `ParamDB::instance().clear()` を呼びテスト間の汚染を防ぐ。

## テスト仕様

`tests/test_params.cpp`:

1. **roundtrip**: `set("B.d0", 0.074)` → `get("B.d0", -1) == 0.074`、
   `all().size() == 1`、`erase` 後 `get("B.d0", -1) == -1`
2. **fallback**: 未設定キー `get("P.d0", 12.5) == 12.5` (bit 一致)
3. **拡散係数の直接検証**: `d = find_dopant("B")`、
   `D_base = dopant_diffusivity(*d, 1273.15, 1.0)` →
   `set("B.d0", 2 * d->d0)` → 新しい D について
   `(D_new - D_base) == d->d0 * exp(-e0/kT)` を相対差 < 1e-12 で確認
   (中性項のみ倍増。B は d0 と dp を持つため全体は 2 倍にならない —
   項単位で検証する)
4. **プロファイル広がりのスケーリング**: 6×6×24 箱メッシュ、B デルタ的
   初期プロファイル (表層 2 セル層に 1e18)。equilibrium `diffuse`
   (T=1000 °C =1273.15 K, time=600 s, field_enh=false) を
   (i) 上書きなし、(ii) `B.d0` と `B.dp` を**両方 2 倍**、で実行し、
   濃度重心まわりの標準偏差 σ = sqrt(Σ C·(z−z̄)²·V / Σ C·V) の比
   `σ_ii / σ_i` が **sqrt(2) に対し相対差 < 5%** (√(2Dt) スケーリング)
5. **未知キーの無害性**: `set("bogus.key", 1.0)` 後に diffuse が
   従来どおり完走し結果が上書きなしと bit 一致

`python/test_comprehensive.py` に追加:

- roundtrip: `sim.set_param("B.d0", 0.05)` → `sim.get_param("B.d0") == 0.05`
- **TED 抑制**: implant(damage=True) → `sim.set_param("B.fi", 0.0)` →
  `sim.diffuse(time=1, temp=950, ted=True)` の B プロファイル σ が、
  同条件で ted=False (equilibrium) の σ に対し**比 < 1.1** (fi=0 で
  TED 増速が消えることの検証)。比較用に SimState を 2 つ (または
  Simulation 2 インスタンス) 使い、テスト末尾で
  `sim.set_param` の上書きを `_c.proc_set_param` ではなく
  新規 `proc_get/erase` — 簡便には各テスト冒頭で不要な上書きが無い前提で
  よい (Python 側に clear は公開しない。上書きした値はテスト内で
  元値に戻す: `sim.set_param("B.fi", 1.0)`)

## 完了条件 (DoD)

- [ ] `ParamDB` (get/set/erase/clear/all, シングルトン, ロックなしの
      制約を文書化) 実装、CMakeLists.txt に src/param_db.cpp 追加
- [ ] 上表の全キーが消費箇所で ParamDB 経由になっている
      (kDopants は無改変、ガード条件は上書き後の値で判定)
- [ ] `proc::set_param/get_param/list_params` + pybind 3 本 +
      `Simulation.set_param/get_param` + Python テスト
- [ ] `tests/test_params.cpp` を foreach に登録、上記 5 テスト全 PASS、
      既存テスト (上書きなしでは挙動不変のはず) 全 PASS
- [ ] コミットメッセージに `P1-10` を含める

## やらないこと

- キーのスキーマ検証・補完・一覧ドキュメント自動生成 (未知キーは黙って保持)
- 材料別 (SiO₂ 中など) のパラメータ名前空間 (P1-9 で拡張)
- ファイル (JSON/TOML) からの一括ロード、デッキコマンド化
- ロックによるスレッド安全化 (solve 中の set 禁止を文書化するのみ)
- range テーブル (Rp/ΔRp) の上書き (implant の rp=/drp= 引数で既に可能)
