# P3-e: 応力依存物理 (拡散・酸化への応力フィードバック)

## 目的

P2-6 (`proc::mechanics`) が解いた応力場 (`st.fields["sxx"].."sxz"`,
dyn/cm²) を、(1) ドーパント拡散係数と (2) 2D 酸化の界面反応速度へ
フィードバックする。結合は P3_overview.md の確定方針どおり:

1. 拡散: `D → D·exp(−p·V_act/kT)`、p = −tr(σ)/3 (静水圧、圧縮で正)。
   V_act は種別 ParamDB キー。dcell 計算ループ内でセル応力を参照する
   だけの**純粋に局所的な陽的結合** (既存 Picard/nni 外側ループの中で
   自然に収束する。陰的結線・Jacobian 変更は一切しない)。
2. 酸化: 界面反応速度 `k_s → k_s·exp(−σ_nn·V_r/kT)` (法線応力による
   減速; bird's beak 先端の自己制限)。P2-4 の `oxidize_2d` の Robin BC
   係数にフック。
3. 活性化 (C_ss) の応力補正は**第一段階では省略** (overview 確定)。

**opt-in 方式**: 新しい `proc::` 関数・引数は追加しない。ParamDB キー
`"stress.couple"` (既定 0 = 完全 off) で有効化する。off のとき、および
on でも応力場が存在しないとき (mechanics 未実行) は、**ゼロコストの
no-op パス**で現行実装とビット一致でなければならない。

**前提**: P2-6、P2-4 (両方マージ済み)。P1-10 ParamDB。

## 現状コード

- `src/process.cpp` `proc::mechanics` (1953 行付近): FEM を解いた要素
  応力を `maxwell_update` で緩和した後、6 成分を
  `st.fields["sxx"/"syy"/"szz"/"sxy"/"syz"/"sxz"]` に **dyn/cm² で**
  格納する (MPa→dyn/cm² 係数 `k = 1e7`)。引張が正の符号規約。
  応力場はセル場なので、`oxidize_2d` の `realize_geometry_col` が
  行う nearest-old 転写 (`for (auto& [sym, conc] : st.fields)`) で
  幾何更新後も自動的に追従する。
- `src/diffusion.cpp` `DiffusionSolver::step_once` (515 行付近):
  Picard 外側ループ内、nni 計算 (553–564 行) の直後にある
  per-cell diffusivity ループ (568–590 行) が dcell を作る:
  ```cpp
  double dv = dopant_diffusivity(dp, T, nni[i]);
  if (o.field_enh && dp.type != DopType::neutral) { ... dv *= ...; }
  dcell[s][i] = dv;
  ```
  重要: この dcell は、その後の **classic 逐次パス (614 行〜) と PA-3
  species-parallel パス (771 行〜) の分岐より前に一度だけ**計算され、
  S-3 Newton パス (`o.use_newton`) も同じ dcell を凍結値として使う。
  したがってこのループに乗算するだけで全 3 パスに同時に効く。
- `src/diffusion.cpp` `step_once_ted` (1043 行〜): TED 版は独自の
  dcell ループ (1208–1232 行、`dopant_diffusivity` → field_enh →
  CI/CV スケール) を持つ。**ここにも同じ因子が必要** (フック 2 箇所)。
- `include/cprocess/diffusion.hpp` `DiffuseOpts` (21 行〜): S-3/S-5/PA
  の opt-in フラグが並ぶ。ここにセル圧力配列ポインタを追加する。
- `src/process.cpp` `proc::diffuse` (1830 行) / `proc::diffuse_ted`
  (1862 行): `DiffusionSolver` を構築し `run`/`run_ted` を呼ぶ唯一の
  結線点。圧力配列はここで st.fields から作る。
- `src/process.cpp` `oxidize_2d` (1521 行〜): N=10 サブステップ。各
  ステップで `ks = BA_cm_s * kNOx / c_gas` (スカラー) を計算し
  `solve_oxidant(mesh, d_cell, active_mask, si_mask, ks, c_gas, log)`
  に渡す。
- `src/oxidant_solver.cpp` `solve_oxidant` (53 行〜): Robin 面は
  「active(酸化膜/窒化膜) セル ↔ silicon セル」間の tet 面 fi ごとに
  `t_robin = 1/(1/t_diff + 1/(ks·area))` を作る (116–137 行、owner 側
  と neigh 側の 2 分岐)。ks はスカラー引数。ここを面別係数に拡張する。
- `include/cprocess/materials.hpp`: `kBoltzmannEv = 8.617333262e-5`
  はあるが erg 単位のボルツマン定数はない (本タスクで追加)。
- P1-10: `proc::set_param`/`get_param` は pybind
  (`python/_cprocess.cpp` 684 行) と `Simulation.set_param`
  (`python/cprocess/simulation.py` 445 行) で既に Python へ露出済み。
  **本タスクの Python 表面は不要** (新 proc:: 関数を作らないため、
  CLAUDE.md の必須ラッパー規則の対象が発生しない)。

## 実装手順

1. **定数とパラメータ**。`include/cprocess/materials.hpp` に
   ```cpp
   constexpr double kBoltzmannErg = 1.380649e-16;  // erg/K
   constexpr double kOmegaSi = 2.0e-23;            // Si 原子体積 [cm^3]
   ```
   を追加。ParamDB キー (すべて `ParamDB::get(key, default)` で参照):
   | キー | 既定値 | 意味 |
   |---|---|---|
   | `stress.couple` | 0.0 | 0=off (完全現行互換)、非 0=on |
   | `stress.vact.B` | 3.4e-24 | B の活性化体積 0.17·Ω_Si [cm³] |
   | `stress.vact.P` | 2.0e-24 | P: 0.10·Ω_Si [cm³] |
   | `stress.vact.As` | 6.0e-24 | As: 0.30·Ω_Si [cm³] |
   | `stress.vact.<他>` | 0.0 | 未定義種は無補正 (因子 exp(0)=1) |
   | `stress.vr` | 1.5e-23 | 酸化反応体積 V_r ≈ 0.75·Ω_Si [cm³] |

   `src/materials.cpp` にヘルパを追加 (B/P/As の既定表 + ParamDB
   オーバーライド、他種は 0):
   ```cpp
   double stress_activation_volume(const Dopant& d);  // [cm^3]
   ```

2. **DiffuseOpts への圧力配列** (`include/cprocess/diffusion.hpp`):
   ```cpp
   // P3-e: per-cell hydrostatic pressure p = -tr(sigma)/3 [dyn/cm^2],
   // size nc. nullptr (default) = coupling off -> dcell loop untouched,
   // bit-identical to pre-P3-e behavior.
   const std::vector<double>* pressure = nullptr;
   ```
   solver 側には SimState を持ち込まない (層分離維持)。

3. **dcell フック (2 箇所)**。`step_once` の dcell ループ (現 577–589
   行) と `step_once_ted` の dcell ループ (現 1217–1230 行) で、Si
   セルの dv 確定直後 (`dcell[s][i] = dv;` の直前、step_once_ted では
   CI/CV スケール適用の**後**) に:
   ```cpp
   if (o.pressure) {  // per-species vact は種ループ先頭で 1 回だけ取得
     const double arg = -(*o.pressure)[i] * vact_s / (kBoltzmannErg * T);
     dv *= std::exp(std::max(-30.0, std::min(30.0, arg)));
   }
   ```
   `vact_s = stress_activation_volume(dp)` は s ループ先頭で hoist。
   `vact_s == 0` の種はループ前に `if (o.pressure && vact_s != 0.0)`
   で分岐ごとスキップ (neutral 種・I/V 等は自動的に無補正)。
   **dcell は classic/PA-3/Newton 分岐より前に共有計算されるため、
   この 2 箇所だけで 3 パスすべてに効くこと**をコードコメントに明記。

4. **proc 層の結線** (`src/process.cpp` `diffuse` と `diffuse_ted`):
   `DiffusionSolver` 構築前に
   ```cpp
   std::vector<double> pressure;  // 関数ローカル、opts はコピーして渡す
   DiffuseOpts o = opts;
   const auto sx = st.fields.find("sxx"), sy = st.fields.find("syy"),
              sz = st.fields.find("szz");
   const std::size_t nc = st.mesh.cells.size();
   if (ParamDB::instance().get("stress.couple", 0.0) != 0.0 &&
       sx != st.fields.end() && sy != st.fields.end() &&
       sz != st.fields.end() && sx->second.size() == nc &&
       sy->second.size() == nc && sz->second.size() == nc) {
     pressure.resize(nc);
     for (std::size_t i = 0; i < nc; ++i)
       pressure[i] = -(sx->second[i] + sy->second[i] + sz->second[i]) / 3.0;
     o.pressure = &pressure;
     if (log) *log << "[diffuse] stress coupling ON (p in "
                   << ... min/max ... << " dyn/cm^2)\n";
   }
   ```
   条件を満たさなければ `o.pressure` は nullptr のまま = 完全 no-op。
   サイズ不一致 (mechanics 後に refine 等でメッシュが変わった) も
   黙って no-op とし、log に 1 行警告を出す。

5. **酸化の Robin 係数フック**。
   (a) `include/cprocess/oxidant_solver.hpp` / `src/oxidant_solver.cpp`:
   `solve_oxidant` の d_cell 版オーバーロードに末尾引数
   `const std::vector<double>* ks_face_scale = nullptr` (サイズ nf、
   面 fi の k_s 乗数) を追加。RFace 生成の 2 分岐 (116–137 行) で
   `const double ks_f = ks * (ks_face_scale ? (*ks_face_scale)[fi] : 1.0);`
   を `t_robin` 式に使う。nullptr 時は完全に現行式。
   (b) `oxidize_2d` (src/process.cpp): サブステップループ内、
   `solve_oxidant` 呼び出し前に、結合 on かつ応力 6 成分場が現メッシュ
   サイズで存在する場合のみ:
   - **サブステップ毎に `proc::mechanics(st, temp_k, dt_step, log)` を
     再実行**して現メッシュ上の応力を更新する (幾何が毎ステップ変わる
     ため。静的場の流用ではなく再解が確定仕様。N=10 回の FEM 解は
     テストメッシュ規模では許容コスト)。結合 off なら一切呼ばない。
   - 全面 fi について Si 側セル c_si (owner/neigh のうち `si_mask` の
     側) の応力テンソルから法線応力を採る:
     ```cpp
     n̂ = f.S / |f.S|;   // 面法線
     σ_nn = n̂ᵀ σ(c_si) n̂
          = nx²·sxx + ny²·syy + nz²·szz
            + 2(nx·ny·sxy + ny·nz·syz + nx·nz·sxz);  // dyn/cm², 引張正
     scale[fi] = exp(min(0, σ_nn) · V_r / kT);       // 圧縮のみ減速、≤ 1
     ```
     指数は −30 でクランプ。Robin 面以外の scale は 1.0 のままで良い
     (solve_oxidant は Robin 分岐でしか読まない)。
   - 引張下の増速はしない (`min(0, σ_nn)`; Kao 型の減速モデルのみ、
     符号規約とともにコードコメントに明記)。

6. **テスト追加** (下記) と `CMakeLists.txt` の
   `foreach(t ... rcm)` に `stress_physics` を追加。

7. Python: 新 proc:: 関数なしのため pybind/`simulation.py` 変更なし。
   `python/test_comprehensive.py` にスモーク 1 本のみ追加:
   `sim.set_param("stress.couple", 1)` → `mechanics` → `diffuse` が
   例外なく走り、`set_param("stress.couple", 0)` の結果と B の spread
   が変化する (方向は圧力符号に依存するので「差が生じる」ことのみ)。
   テスト末尾で `stress.couple` を 0 に戻す。

## テスト仕様

`tests/test_stress_physics.cpp` 新設。メッシュは既存テスト流の
`mesh_box(st, 0, 0.2e-4, 0, 0.2e-4, 0, 0.5e-4, 4, 4, 50)` を基本とする。

1. **解析検証 (一様静水圧)**: B Gaussian (rp=0.1 µm, 1e18) を埋込み、
   `st.fields["sxx"]=["syy"]=["szz"] = −1e10` (一様 1 GPa 圧縮)、
   `sxy=syz=sxz=0` を**手で**セット (mechanics は呼ばない)。
   `stress.couple=1` で `diffuse(t, 1000 °C)`。
   f = exp(−p·V_act(B)/kT)、p = 1e10 dyn/cm²、T = 1273.15 K
   (f ≈ 0.824)。比較対象: 応力なし・時間 `t·f` の diffuse。
   両者の深さプロファイルの spread (test_ted.cpp の `profile_spread`
   流用) の相対差 **< 2%** (√(D t) 等価性による解析チェック)。
2. **応力場なしでビット一致**: `stress.couple=1` だが応力場を
   持たない状態で diffuse / diffuse_ted → 結合なし
   (`stress.couple=0`) の同一シーケンスと B 濃度ベクトルが
   **全セルで完全一致 (==、許容 0)**。no-op パスは同一コード経路の
   ため丸め差も生じないこと (overview の < 1e-12 より強い基準を採用)。
3. **couple=0 でビット一致**: 応力場を 1 GPa でセットしても
   `stress.couple=0` (既定) なら現行結果と完全一致 (==)。
4. **PA-3 / Newton パス整合**: テスト 1 の設定 + 3 種 (B/P/As) で
   `species_parallel=1` と classic (`species_parallel=-1`) の結果が
   相対差 < 1e-12、`use_newton=true` でも spread がテスト 1 の解析値と
   < 2%。
5. **LOCOS bird's beak 短縮 (定性)**: test_ox2d.cpp の LOCOS 設定
   (窒化膜マスク付き wet 酸化) を流用。A: `stress.couple=0`、
   B: `stress.couple=1` (oxidize_2d 内部で mechanics が回る)。
   マスク下 taper 長 (test_ox2d と同じ計測) が **B で A の 70–95%**
   (5–30% 短縮)。open-field 膜厚の変化は < 10% (減速は圧縮の強い
   先端に局在すること)。
6. 既存テスト全 PASS。特に **test_diffusion / test_ted /
   test_species_parallel / test_newton / test_adaptive_dt / test_ox2d /
   test_fem / test_oed** は本タスクの触る経路の回帰として必須。

Python: `python/test_comprehensive.py` に手順 7 のスモーク。
`python3 python/test_comprehensive.py` 全 PASS。

## 完了条件 (DoD)

- [ ] `stress.couple` 既定 0 で全経路が現行とビット一致 (テスト 2, 3)
- [ ] dcell フックが step_once / step_once_ted の 2 箇所に入り、
      classic / PA-3 / Newton の 3 パスで同一に効く (テスト 4)
- [ ] 一様 1 GPa で B の D スケーリングが解析値と < 2% (テスト 1)
- [ ] `solve_oxidant` の `ks_face_scale` が nullptr 時に現行と同一、
      LOCOS taper 5–30% 短縮 (テスト 5)
- [ ] `tests/test_stress_physics.cpp` + CMakeLists 登録、
      C++ 全テスト PASS、Python テスト PASS
- [ ] 新 proc:: 関数ゼロ (ParamDB のみ)。コミットメッセージに `P3-e`

## やらないこと

- 活性化 (C_ss) の応力補正 (overview で第一段階省略と確定)
- 完全異方性テンソル拡散 (D をテンソル化しない — 静水圧スカラー補正のみ)
- 応力依存の偏析係数・界面輸送 (P1-4 の m/h は不変)
- 拡散→応力の逆結合 (組成起因の固有ひずみは P3-f)、
  拡散ステップ内での mechanics 自動再実行 (ユーザが明示的に呼ぶ;
  自動再解は oxidize_2d のサブステップ内のみ)
- 点欠陥 (I/V/C311) 自身の輸送への応力補正 (vact 既定 0 のまま)
- 1D `oxidize()` (Deal-Grove 解析式) への応力フック — 2D 専用
- 新しい proc:: 関数 / デッキコマンド / Simulation メソッド
