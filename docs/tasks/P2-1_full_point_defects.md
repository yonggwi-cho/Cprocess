# P2-1: フル点欠陥モデル (I + V 連立 + {311} クラスタ)

## 目的

現行の TED (`DiffusionSolver::run_ted`, src/diffusion.cpp) は自己格子間原子
(I) の**過剰分 ψ のみ**を単一の線形場として解き、{311} クラスタは
`kSmax = 3.0e3` の S 上限で近似している。本タスクで置き換える:

- 格子間原子 **C_I** と空孔 **C_V** の 2 場 (**絶対濃度**、過剰分ではない) を
  I-V バルク再結合つきで連立輸送する
- 不動の **{311} クラスタ場 C_311** を明示的に導入し、I の捕獲/放出で
  TED の持続時間 (リザーバ効果) を物理的に再現する
- ドーパント増速を `(1−fi) + fi·(C_I/C_I*) + fv·(C_V/C_V*)` (fv = 1−fi) に
  拡張し、Sb 等の空孔支配種にも正しい応答を与える

**P1-2 (MC 損傷 → I シード) と P1-10 (ParamDB) の完了が前提。**
非線形結合はオペレータ分割 (反応サブサイクル) で扱う。S-2/S-3 が完了すれば
将来ブロック陰解法へ移行できるが、本タスクは分割法のみとする。

## 現状コード

- `src/diffusion.cpp` — `run_ted(fields, psi, bcs, opts)`:
  ステップごとに (1) ψ を陰解法 1 ステップ (拡散 + 一次反応 `k_rec`、
  表面 zmax で ψ=0 Dirichlet)、(2) `S = 1 + min(psi/cstar, kSmax)` で
  ドーパント D を `(1−fi) + fi·S` 倍して Picard で進める。
  `assemble(dcell, cold, bcface, cgrad, dt, reaction, nonortho, rhs, grad)` は
  一次反応係数 `reaction` [1/s] を対角に加算できる (ψ の陰解法で使用中)。
- `include/cprocess/diffusion.hpp` — `run_ted` 宣言、`DiffuseOpts`。
- `src/materials.cpp` — `interstitial_cstar / interstitial_diffusivity /
  interstitial_recomb_rate` (Arrhenius)、`Dopant::fi`。
- `src/process.cpp` — `seed_interstitials` ("+1"): 注入増分を `st.fields["I"]`
  に**過剰濃度として**加算。`proc::diffuse_ted` が `st.fields["I"]` を ψ として
  `run_ted` に渡す。P1-2 完了後は MC の Kinchin-Pease 損傷が同じ "I" 場
  (過剰濃度) をシードする。
- `tests/test_ted.cpp` — 3 テスト (増速 >1.3x、transient 減衰、B>As)。
- P1-10 の ParamDB: `sim.set_param("key", value)` / C++ 側
  `params.get("key", default)` で物性値を上書きできる機構 (P1-10 仕様参照)。

## 実装手順

1. **パラメータ (materials.cpp + ParamDB)**。以下を ParamDB キーとして
   登録し、既定値は教科書オーダーの値をハードコードする
   (すべて `params.get(key, default)` 経由で読む):

   | キー | 意味 | 既定値 |
   |------|------|--------|
   | `pd.ci_star.pre` / `pd.ci_star.e` | C_I* = pre·exp(−E/kT) | 3.0e27 cm⁻³ / 3.7 eV (現行 `interstitial_cstar`) |
   | `pd.cv_star.pre` / `pd.cv_star.e` | C_V* = pre·exp(−E/kT) | 1.0e27 cm⁻³ / 3.6 eV |
   | `pd.di.pre` / `pd.di.e` | D_I | 5.0e-2 cm²/s / 1.77 eV (現行) |
   | `pd.dv.pre` / `pd.dv.e` | D_V | 1.0e-3 cm²/s / 1.8 eV |
   | `pd.kbulk.factor` | k_bulk = factor·4π·a_Si·(D_I+D_V)、a_Si=2.35e-8 cm | 1.0 |
   | `pd.c311.ktrap.factor` | k_trap = factor·4π·a_Si·D_I | 1.0 |
   | `pd.c311.nu0` | 放出試行頻度 ν0 | 1e13 s⁻¹ |
   | `pd.c311.eb` | {311} 結合エネルギー E_b | 3.6 eV |
   | `pd.damage.v_fraction` | 初期 V/I シード比 | 0.9 |

   `materials.hpp` に関数を追加 (ParamDB を引数で受ける):
   ```cpp
   struct PointDefectParams {
     double ci_star, cv_star;   // cm^-3 (at T)
     double d_i, d_v;           // cm^2/s
     double k_bulk;             // cm^3/s  (bimolecular I-V recombination)
     double k_trap;             // cm^3/s  (I capture by {311})
     double k_emit;             // 1/s     = nu0 * exp(-eb/kT)
   };
   PointDefectParams point_defect_params(double temp_k, const ParamDB& db);
   ```
   既存の `interstitial_cstar/diffusivity/recomb_rate` は残す (他の呼び出し元
   互換のため) が、run_ted 内では使わなくなる。

2. **支配方程式** (Si セルのみ、非 Si は凍結 = 現行 mask 方式):
   ```
   ∂C_I/∂t = ∇·(D_I ∇C_I) − k_bulk (C_I C_V − C_I* C_V*)
             − k_trap C_I · H(C_I − C_I*) · (C_I − C_I*) / C_I  + k_emit C_311
   ∂C_V/∂t = ∇·(D_V ∇C_V) − k_bulk (C_I C_V − C_I* C_V*)
   dC_311/dt = k_trap · max(C_I − C_I*, 0) − k_emit · C_311      (不動、拡散なし)
   ```
   つまり {311} 捕獲項は具体的に `k_trap · max(C_I − C_I*, 0)`
   [cm⁻³/s、k_trap は cm³/s なので max(...) に C_311 相当の濃度次元を
   持たせないよう **k_trap·max(C_I − C_I*,0) をそのまま体積レートとする —
   実装では rate_trap[i] = k_trap * std::max(CI[i] - ci_star, 0.0) *
   std::max(CI[i], 0.0) / std::max(CI[i], 1.0)** ではなく単純化して
   `rate_trap = k_trap_eff * max(C_I − C_I*, 0)` (k_trap_eff [1/s] =
   k_trap · C_ref, C_ref = 1e22 cm⁻³ 固定) とする。放出は
   `rate_emit = k_emit * C_311`。C_I への正味項は `−rate_trap + rate_emit`、
   C_311 へは `+rate_trap − rate_emit`。**この形をコードコメントに明記する。**

3. **`run_ted` の内部を置換** (シグネチャは変えない):
   ```cpp
   void run_ted(std::vector<SpeciesField>& fields, std::vector<double>& psi,
                const std::vector<DirichletBC>& bcs, const DiffuseOpts& o);
   ```
   - (a) **初期化**: 入力 `psi` は過剰 I 濃度 (P1-2/+1 シード)。内部で
     `CI[i] = ci_star + psi[i]`、`CV[i] = cv_star + v_fraction * psi[i]`
     (`pd.damage.v_fraction`=0.9 — バルク再結合が大半を食う)、
     `C311[i] = 0` を組み立てる。C_311 は `run_ted` 呼び出し間で保持する
     必要があるため、**proc 層で `st.fields["C311"]` を渡せるよう**
     `run_ted` に第 3 の場を追加したオーバーロードを新設:
     ```cpp
     void run_ted(std::vector<SpeciesField>& fields, std::vector<double>& psi,
                  std::vector<double>& c311, const std::vector<DirichletBC>& bcs,
                  const DiffuseOpts& o);
     ```
     旧シグネチャはローカル c311 (ゼロ初期化、破棄) で新版を呼ぶ薄い転送に
     する。`proc::diffuse_ted` は `st.fields["C311"]` (無ければ作成) を渡す。
   - (b) **境界条件**: C_I = C_I*、C_V = C_V* を zmax パッチで Dirichlet
     (現行 ψ=0 と等価な絶対濃度版)。他は零フラックス。
   - (c) **タイムステップのオペレータ分割** (各グローバル dt につき):
     1. **拡散 (陰)**: C_I を `assemble(dI, CI_old, ci_bcface, CI, dt, 0, ...)`
        で 1 ステップ、C_V も同様 (反応項はここに入れない)。
     2. **反応 (陽・サブサイクル)**: セルごとに独立の 3 変数 ODE を
        前進オイラーでサブステップ。刻みは
        ```cpp
        const double k_rate = k_bulk * std::max(CI[i]*1.0, CV[i]*1.0); // 1/s の目安
        // 全セル一括: dt_react = min(dt, 0.1 / max_i(k_bulk*max(CI,CV) + k_trap_eff + k_emit))
        ```
        具体的に: `rate_max = max_i( k_bulk*max(CI[i],CV[i]) + k_trap_eff + k_emit )`、
        `dt_react = std::min(dt, 0.1 / std::max(rate_max, 1e-30))`、
        `nsub = ceil(dt / dt_react)` で `dt/nsub` 刻み。各サブステップで
        ```
        R  = k_bulk (CI·CV − ci_star·cv_star)
        Tr = k_trap_eff · max(CI − ci_star, 0)
        Em = k_emit · C311
        CI   += dts (−R − Tr + Em);  CV += dts (−R);  C311 += dts (Tr − Em)
        ```
        後に負値クランプ (max(v,0))。
     3. **ドーパント (陰、Picard)**: 現行の Picard ループを流用し、増速のみ
        ```cpp
        const double fv = 1.0 - dp.fi;
        dv *= (1.0 - dp.fi) + dp.fi * (CI[i] / ci_star) + fv * (CV[i] / cv_star);
        ```
        …ではなく **(1−fi) 項は削除**し、正しくは
        `dv *= dp.fi * (CI[i]/ci_star) + (1.0 - dp.fi) * (CV[i]/cv_star)`
        とする (平衡 CI=CI*, CV=CV* でちょうど 1 になる)。S 上限
        `kSmax` は {311} が物理的に緩衝するため**撤廃**するが、数値安全の
        ため合計スケールを 1e4 でクランプ (`std::min(scale, 1e4)`)。
   - (d) 終了時 `psi[i] = CI[i] − ci_star` (負も許す→0 クランプ) に書き戻し、
     `st.fields["I"]` の意味 (過剰濃度) と proc:: API の互換を保つ。
     `st.fields["V"]` にも `CV − cv_star` を書き出す (proc 層で作成)。
   - (e) ログ: 各報告行に `Smax=CImax/CI*`、`C311max` を出す。
4. **proc:: / Python**: `proc::diffuse_ted` のシグネチャは不変なので
   新規バインディングは不要。ただし `st.fields` に "V"/"C311" が増えるため
   Python の `sim.field("V")`, `sim.field("C311")` が動くことをテストで確認。
   ParamDB 経由の調整 (`sim.set_param("pd.c311.eb", 3.8)`) の smoke テストを
   追加する。
5. `interstitial_recomb_rate` を使っていた箇所 (run_ted のみ) を整理。

## テスト仕様

`tests/test_ted.cpp` を改修 + `tests/test_point_defects.cpp` (新規、
CMakeLists の foreach に `point_defects` 追加)。

1. **平衡安定性** (test_point_defects): 損傷なし (`damage=false`) の B 注入
   → `diffuse_ted(30 min, 900 °C)`。終了後 `max_i |CI[i] − ci_star| / ci_star
   < 0.01` (psi 経由: `max(fields["I"]) < 0.01*ci_star`)。同条件の
   `proc::diffuse` (plain) と B の spread が相対差 < 2% で一致。
2. **TED 過渡の再現** (test_ted 改修): 既存テスト 1〜3 は維持。合格基準で
   正当に変わりうるもの: (a) 増速倍率の絶対値 (kSmax 撤廃 + {311} 緩衝で
   変化しうる — `spread_ted > 1.3 * spread_eq` は維持すること)、
   (b) テスト 2 の d2 < d1 は維持、(c) B > As は維持。数値ログの文言
   (`Smax=` の定義) は変わってよい。
3. **{311} リザーバ効果** (test_point_defects): 同一の損傷つき B 注入 state
   を 2 つ用意し、片方は `set_param("pd.c311.ktrap.factor", 0.0)` で
   クラスタ無効。900 °C で t=30 s と t=300 s 時点の spread を比較:
   - クラスタ有効: `spread(300s) − spread(30s) > 0.2 * (spread(30s) − spread(0))`
     (後期にも増速が持続)
   - クラスタ無効: `spread(300s) − spread(30s) < spread(30s) − spread(0)`
     (増速が前半に集中して早く死ぬ)
   - かつ有効時の後期増分 > 無効時の後期増分。
4. **I-V 再結合** (test_point_defects): メッシュ中央帯に
   `st.fields["I"] = st.fields["V"] = 1e19` (psi/v 過剰を直接 set_field)、
   ドーパントなしで `diffuse_ted(60 s, 1000 °C)` → 終了時
   `max(fields["I"]) < 0.1 * 1e19` かつ `max(fields["V"]) < 0.1 * 1e19`
   (両者とも平衡へ減衰)。
5. **反応サブサイクルの妥当性**: テスト 4 を `o.dt` を 2 倍にして再実行し
   spread/残留値が相対差 < 5% (dt_react 制御が効いている)。

Python (`python/test_comprehensive.py`): 損傷つき implant → `diffuse(ted=True)`
→ `sim.field("V")` / `sim.field("C311")` が取得でき、`C311` の総量が非負。
`set_param("pd.c311.nu0", 1e14)` 後の再実行が完走する。

## 完了条件 (DoD)

- [ ] `run_ted` が C_I/C_V/C_311 の分割解法に置換され、旧シグネチャ +
      `proc::diffuse_ted` が互換動作 (psi = C_I − C_I*)
- [ ] 全パラメータが ParamDB キー経由 (上表のキー名どおり)
- [ ] 上記 C++ テスト 5 本 + 改修 test_ted が PASS、CMakeLists 登録
- [ ] Python テスト追加、`python3 python/test_comprehensive.py` PASS
- [ ] `cmake --build build -j$(nproc) && ctest --test-dir build` 全 PASS
- [ ] コミットメッセージに `P2-1` を含める

## やらないこと

- 5-stream (dopant-I / dopant-V ペア場の明示輸送) — 将来 S-2/S-3 完了後
- ブロック陰解法・Newton-Krylov (S-2/S-3)。本タスクは分割 + 陽反応のみ
- BIC / As-V ドーパントクラスタ (P2-2)
- OED の酸化界面 I 注入 (P2-3)
- {311} のサイズ分布 / Ostwald ripening (単一有効場のみ)
- 転位ループ、荷電点欠陥 (I⁰/I⁺/I⁻ 区別なし)
