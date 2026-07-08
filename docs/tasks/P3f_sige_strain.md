# P3-f: SiGe / 歪みエンジニアリング (Ge 組成 → 固有ひずみ・ni 補正)

## 目的

P2-8 で導入済みの Ge 場 (`st.fields["Ge"]`, cm⁻³, DopType::neutral の
不動マーカー種) から物理を生成する。P3_overview.md の確定方針:

1. **格子失配固有ひずみ** (Vegard 則): ε0 = 0.042·x_Ge、
   x_Ge = C_Ge / 5e22 (`kNSi`、materials.hpp に P2-3 で公開済み)。
   これを P2-6 `proc::mechanics` の固有ひずみ荷重に、既存の
   熱失配・真性膜応力項と**並列に加算**する (セル毎・組成由来項)。
2. **バンドギャップ縮小**: ΔEg = 0.4·x_Ge [eV] →
   ni_SiGe = ni·exp(ΔEg/2kT)。拡散の電界増速・濃度依存 D が使う
   n/ni 評価点 (nni 計算) をセル毎 ni に置換する。
3. **拡散への歪み補正は静水圧成分のみ、P3-e の機構経由**。P3-f は
   ひずみ源 (ε0) を供給するだけで、輸送補正そのものは
   「mechanics で応力を解く → `stress.couple` + `DiffuseOpts::pressure`
   (P3e_stress_physics.md 手順 2–4)」が行う。本タスクで拡散係数式には
   一切触らない。

**前提**: P2-6、P2-8、**P3-e (先に実装必須** — 本仕様は P3-e が追加する
`DiffuseOpts::pressure` / `stress.couple` 機構の存在を前提とする)。

## 現状コード

- `src/materials.cpp` (107–125 行): Ge エントリ。DopType::neutral、
  d0=0 (全材料で D=0、不動)。**neutral のため diffusion.cpp の nnet
  集計 (`if (type == DopType::neutral) continue;`) から除外済み** —
  Ge は nnet 自体を乱さず、ni の値だけを変えるのが本タスクの効き方。
- `include/cprocess/materials.hpp`: `kNSi = 5.0e22` (P2-3)、
  `kBoltzmannEv = 8.617333262e-5`。
- `src/materials.cpp` `ni_si(temp_k)` (148 行): 温度のみの関数。
  シグネチャは変えない (セル毎補正は diffusion.cpp 側で乗算する)。
- `src/diffusion.cpp` — ni 評価点は **2 箇所**:
  - `step_once` (515 行〜): `const double ni = ni_si(T);` (527 行) を
    Picard ループ内の nni 計算 (553–564 行) が
    `cc = nnet/(2·ni); nni[i] = cc + sqrt(cc²+1)` として使う。
  - `step_once_ted` (1043 行〜): 同型の nni ループ (1196–1207 行)。
  nni[i] はその直後の dcell ループ (`dopant_diffusivity(dp, T, nni[i])`
  と field_enh 因子) に入る。dcell は classic / PA-3 species-parallel /
  S-3 Newton の分岐**前**に共有計算されるため (P3-e 仕様の現状コード欄
  参照)、nni ループ 2 箇所の修正で全 3 パスに効く。
- `src/process.cpp` `proc::mechanics` (1953 行〜): 固有ひずみ構築
  ループ (1973–1986 行) が
  `eps0[c*6+k] = e_thermal + e_intrinsic (k=0..2)` を書く。ここが
  組成由来 ε0 の加算点。`fem_assemble(st.mesh, E_cell, nu_cell, eps0)`
  (src/fem.cpp 62 行) は eps0_cell を要素荷重 `f_e += V·Bᵀ·D·ε0`
  (141–150 行) と要素応力 `σ = D(Bu − ε0)` (fem_element_stress,
  223–228 行) の両方に既に配管済み — **fem.cpp は無変更**。
- `st.fields["Ge"]` は implant/init で作られるセル場。mechanics は
  SimState を持つので直接参照できる。DiffusionSolver は SimState を
  持たないが、Ge は `fields` (SpeciesField 列) に dopant として入って
  くるので solver 内でシンボル "Ge" を探せば良い。
- P1-10 ParamDB: Python へは `set_param`/`get_param` で露出済み
  (P3-e 仕様と同じ理由で本タスクも **Python 表面は不要**)。

## 実装手順

1. **ParamDB キー** (全て `ParamDB::get`):
   | キー | 既定値 | 意味 |
   |---|---|---|
   | `sige.couple` | 1.0 | 0 で SiGe 物理を完全 off (kill switch) |
   | `sige.eps0_coef` | 0.042 | Vegard 係数: ε0 = coef·x_Ge |
   | `sige.dEg_coef` | 0.4 | ΔEg = coef·x_Ge [eV] |

   既定 on で良い理由: Ge 場が無い/ゼロなら x_Ge = 0 →
   `eps0 += 0.042·0 = 0`、`exp(0.4·0/2kT) = 1.0` で**厳密に**現行と
   一致する (overview の「Ge=0 で全結果が現行と一致」を恒等式で満たす)。
   Ge 場が存在しないときは参照ループ自体をスキップし、ゼロコスト。

2. **固有ひずみ** (`src/process.cpp` `proc::mechanics`)。
   eps0 構築ループ (1973–1986 行) の後に:
   ```cpp
   const auto ge = st.fields.find("Ge");
   const double eps_coef = db.get("sige.eps0_coef", 0.042);
   if (db.get("sige.couple", 1.0) != 0.0 && ge != st.fields.end() &&
       ge->second.size() == static_cast<std::size_t>(nc)) {
     for (int c = 0; c < nc; ++c) {
       if (!is_silicon(mat_name[c]) && mat_name[c] != "poly" &&
           mat_name[c] != "polysilicon") continue;  // Ge は Si 格子内のみ
       const double xge = std::clamp(ge->second[c] / kNSi, 0.0, 1.0);
       const double e_misfit = eps_coef * xge;  // 等方体積ひずみ近似
       for (int k = 0; k < 3; ++k) eps0[c * 6 + k] += e_misfit;
     }
     if (log) *log << "mechanics: SiGe eigenstrain ON (max x_Ge=...)\n";
   }
   ```
   熱項と同じ「等方 (静水) 固有ひずみ」形式を採る (P2-6 の真性膜応力と
   同じ文書化済み簡略化 — 実際のミスフィットは面内 2 軸だが、
   ε0 構築の一貫性と P3-e の静水圧結合 (tr(σ)/3) には等方形で十分)。
   剪断成分 (k=3..5) は 0 のまま。

3. **ni 補正** (`src/diffusion.cpp`、2 箇所)。
   `step_once` と `step_once_ted` の nni ループ:
   - ループ前に Ge の species index を 1 回だけ探す:
     ```cpp
     int s_ge = -1;
     for (int s = 0; s < ns; ++s)
       if (fields[s].dopant->symbol == "Ge") { s_ge = s; break; }
     const double dEg_coef = ParamDB::instance().get("sige.dEg_coef", 0.4);
     const bool sige_on =
         s_ge >= 0 && ParamDB::instance().get("sige.couple", 1.0) != 0.0;
     ```
   - nni ループ内 (Si セル分岐) で:
     ```cpp
     double ni_c = ni;                        // = ni_si(T)
     if (sige_on) {
       const double xge = std::min((*fields[s_ge].conc)[i] / kNSi, 1.0);
       if (xge > 0.0)
         ni_c *= std::exp(dEg_coef * xge / (2.0 * kBoltzmannEv * T));
     }
     const double cc = nnet / (2.0 * ni_c);
     nni[i] = cc + std::sqrt(cc * cc + 1.0);
     ```
   Ge 自身は neutral なので nnet 集計は無変更 (P2-8 の
   `continue` がそのまま効く)。nni[i] が変わることで
   `dopant_diffusivity` の nni 冪項と field_enh 因子の両方に補正が
   伝播する。dcell は分岐前共有 (現状コード欄) なので **classic /
   PA-3 / Newton の全パスで同一に効く**。ni 増加 → 同じ nnet に対し
   |cc| 減少 → 外因性領域で濃度依存 D と電界増速が弱まる、が期待符号。
   Picard 収束判定・ループ構造は無変更 (nni は元々毎パス再計算)。

4. **拡散への歪み補正の結線はしない**。ユーザフロー
   `implant Ge → mechanics(T, dt) → set_param("stress.couple",1) →
   diffuse` で、手順 2 の ε0 が応力場 sxx..sxz を生み、P3-e の
   `DiffuseOpts::pressure` 経由で D·exp(−p·V_act/kT) が掛かる。
   このフローを process.hpp の mechanics / diffuse の API コメントに
   1 行ずつ追記するのみ。

5. **テスト追加** (下記) と `CMakeLists.txt` の `foreach(t ...)` に
   `sige` を追加。

6. Python: 新 proc:: 関数なし → pybind / `simulation.py` 変更なし
   (CLAUDE.md 必須ラッパー規則の対象なし)。
   `python/test_comprehensive.py` にスモーク 1 本: Ge 1e22 を含む
   `mechanics` で sxx が finite かつ SiGe 層で圧縮 (負)、
   `sige.couple=0` に戻して終了。

## テスト仕様

`tests/test_sige.cpp` 新設。メッシュ:
`mesh_box(st, 0, 0.2e-4, 0, 0.2e-4, 0, 0.5e-4, 4, 4, 50)`。

1. **ミスフィットひずみの定量 (自由膨張)**: 全セルに一様
   Ge = 1e22 cm⁻³ (x_Ge = 0.2)、`mechanics(300 K 相当 dT=0, dt=0)`。
   min 面ローラー BC (P2-6) の下で一様等方 ε0 は自由膨張するので、
   FEM 解の最上面ノードの uz / Lz が **0.0084 (= 0.042·0.2) の
   ±2% (相対)** に一致し、von Mises 応力 < 1 MPa (P2-6 テスト 2 と
   同じ判定方式)。※ mechanics が u を返さない場合は、代替として
   「一様 Ge → von Mises < 1 MPa」+「下半分のみ Ge の 2 層構成 →
   SiGe 層の面内応力 σ_xx が圧縮 (負) で
   |σ_xx| = E/(1−ν)·0.0084·(2/3±30%) のオーダー」で判定して良いが、
   0.84%±2% の一次判定は ε0 加算値そのものを
   `fem_element_stress` 入力前に assert する単体チェック
   (eps0[c*6+0] == 0.042·x_Ge を機械精度で) で必ず担保する。
2. **ni 増加と拡散変化の符号 (定性)**: B Gaussian (1e19, rp=0.1 µm)
   を全域一様 Ge=1e22 の有無で `diffuse(30 min, 1000 °C)` 比較。
   (a) Ge あり側で nni 由来の実効 D が下がる方向 (p 型・外因性で
   ni 増 → |cc| 減 → 濃度依存項/電界項の増速が縮む) となり、
   B の spread が Ge なしより**小さい**こと (符号確認のみ、量は不問)。
   (b) 低濃度 (1e15, 真性領域) では両者の spread 相対差 < 1%
   (ni 補正は nnet≪ni では効かない)。
3. **Ge=0 でビット一致**: Ge 場が存在しない状態、および
   `st.fields["Ge"]` が全ゼロの状態の両方で、diffuse / diffuse_ted /
   mechanics の結果 (B 濃度・sxx) が本タスク前と**完全一致 (==)**。
   `sige.couple=0` で Ge=1e22 でも現行と完全一致。
4. **P3-e 連成 (統合・定性)**: 下半分 SiGe (Ge=1e22)・上半分 Si に
   B 一様 1e18 → `mechanics(1000 °C, 0)` →
   `set_param("stress.couple", 1)` → `diffuse(30 min, 1000 °C)`。
   例外なく完走し、`stress.couple=0` の同一フローと B プロファイルに
   有意差 (max 相対差 > 1e-6) が出ること (P3-e が SiGe 起因応力を
   拾っている証跡)。
5. **PA-3 整合**: テスト 2(a) を B/P/As + Ge の 4 種で
   `species_parallel=1` と classic で実行し相対差 < 1e-12。
6. 既存テスト全 PASS。特に **test_fem / test_mechanics /
   test_new_dopants (Ge の不動性) / test_diffusion / test_ted /
   test_species_parallel / test_stress_physics (P3-e)** は回帰必須。

Python: 手順 6 のスモーク、`python3 python/test_comprehensive.py`
全 PASS。

## 完了条件 (DoD)

- [ ] ε0 加算が proc::mechanics の eps0 ループ直後に入り、
      x_Ge=0.2 で 0.0084±2% (テスト 1)
- [ ] ni 補正が step_once / step_once_ted の nni ループ 2 箇所に入り、
      classic / PA-3 / Newton で同一 (テスト 2, 5)
- [ ] Ge=0 / `sige.couple=0` で全結果が現行とビット一致 (テスト 3)
- [ ] P3-e 経由の歪み→拡散フローが端から端まで通る (テスト 4)
- [ ] `tests/test_sige.cpp` + CMakeLists 登録、C++ / Python 全テスト
      PASS。新 proc:: 関数ゼロ。コミットメッセージに `P3-f`

## やらないこと

- 転位・塑性緩和・臨界膜厚 (Matthews–Blakeslee) 判定 — ε0 は常に
  完全整合 (fully strained) を仮定
- Ge 依存の弾性定数 (E, ν は Si のまま; 差は ~10% で第一段階外)、
  Ge 依存の熱膨張係数
- ひずみテンソル異方成分の拡散補正 (静水圧のみ、P3-e 確定方針) と
  バンド分裂による移動度/種別異方性
- Ge 自身の拡散・偏析 (d0=0 の不動マーカーのまま; Si-Ge 相互拡散なし)
- ni_si() のシグネチャ変更や材料表への SiGe 材料 (MatId) 追加 —
  組成はあくまで Si セル上の連続場
- エピタキシャル SiGe 成長コマンド (P3-b で扱う)
- 新しい proc:: 関数 / デッキコマンド / Simulation メソッド
