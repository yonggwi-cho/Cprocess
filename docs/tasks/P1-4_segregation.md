# P1-4: Si/SiO₂ 界面偏析 + dose loss

**依存: P1-6 (酸化統合)。** oxide 領域がメッシュ上に存在すること
(`proc::oxidize` または `proc::deposit("oxide", ...)`) が前提。

## 目的

拡散ソルバーは現在、非 Si セルを凍結し材料界面を零フラックス壁として扱う
(`kBoundOwner`/`kBoundNeigh`)。実際の Si/SiO₂ 界面ではドーパントが偏析
係数 m に従って両相に分配され、特に B は酸化膜側へ吸い出されて Si 中の
dose が減る (dose loss)。本タスクで拡散ソルバーに **3 相偏析の界面
フラックス項**を陰的に組み込み、酸化膜セルを (小さな一定拡散係数で)
解に参加させ、B の dose loss / P・As・Sb の Si 残留を再現する。

## 現状コード

- `include/cprocess/diffusion.hpp` — `DiffusionSolver(const Mesh&,
  std::vector<char> solve_mask, std::ostream* log)`。`FaceKind` は
  `kInactive/kInternal/kBoundOwner/kBoundNeigh`。`FGeom` に二点フラックス幾何。
- `src/diffusion.cpp` — `build()` が面種別・CSR パターン
  (kInternal 面のみ隣接対を追加)・`fslot_`・面彩色 (kInternal は owner+neigh、
  kBoundOwner は owner のみに書く前提の `writes_cells`) を構築。
  `assemble()` が面ループで対角/非対角/rhs を組む。`run()` は種ごとに
  `cg_ilu0` → 失敗時 `bicgstab_ilu0` フォールバック。マスク 0 セルは
  恒等行 (`A=1, rhs=c`) で凍結。
- `include/cprocess/materials.hpp` / `src/materials.cpp` — `Dopant` 構造体
  (Fair モデル係数、`fi`、range テーブル)。`kDopants` テーブル (B/P/As/Sb)。
- `src/process.cpp` — `proc::diffuse` が `silicon_mask(st)` を渡して
  solver を構築。`st.region_material` (タグ→材料名)。

## 実装手順

1. **`Dopant` 構造体の拡張** (`include/cprocess/materials.hpp`):
   ```cpp
   // ── Si/SiO2 interface segregation (P1-4) ──
   // Equilibrium segregation coefficient m(T) = C_si / C_ox at the interface:
   //   m(T) = seg_m0 * exp(-seg_e / kT)
   // Interface transport coefficient (cm/s):
   //   h(T) = seg_h0 * exp(-seg_he / kT)
   double seg_m0 = 10.0, seg_e = 0.0;    // default: Si-favoring, T-independent
   double seg_h0 = 1.0e5, seg_he = 2.0;  // ~1.2e-3 cm/s at 1000 C
   // Diffusivity in SiO2 (plain Arrhenius, cm^2/s): D_ox = dox0*exp(-eox/kT)
   double dox0 = 0.0, eox = 0.0;         // 0 => immobile inside oxide
   ```
   **符号規約を固定する**: m ≡ 平衡時の C_si/C_ox (教科書の偏析係数)。
   m < 1 → 酸化膜側に富む (B)、m > 1 → Si 側に残る (P/As/Sb)。
   界面フラックス (単位面積あたり、**酸化膜 → Si 向き正**):
   ```
   F = h(T) * ( m(T) * C_ox − C_si )      [cm^-2 s^-1]
   ```
   平衡で C_si = m·C_ox。
2. **`kDopants` の既定値** (`src/materials.cpp`、教科書レベルの既定・較正可
   とコメント):
   | 種 | seg_m0 | seg_e (eV) | seg_h0 (cm/s) | seg_he (eV) | dox0 (cm²/s) | eox (eV) |
   |----|--------|-----------|---------------|-------------|--------------|----------|
   | B  | 6.0    | 0.33      | 1.0e5 | 2.0 | 1.23e-4 | 3.39 |
   | P  | 10.0   | 0.0       | 1.0e5 | 2.0 | 0.19    | 4.03 |
   | As | 10.0   | 0.0       | 1.0e5 | 2.0 | 3.7e-2  | 3.70 |
   | Sb | 10.0   | 0.0       | 1.0e5 | 2.0 | 2.6e-2  | 4.00 |
   参考値: B の m(1273 K) = 6.0·exp(−0.33/0.1097) ≈ **0.30**、
   m(1373 K) ≈ 0.37。h(1273 K) ≈ 1.2e-3 cm/s (セル寸法 ~0.01 µm に対し
   界面平衡到達 τ ~ ms、アニール時間より十分速い)。B の
   D_ox(1273 K) ≈ 4.7e-18 cm²/s (30 min で拡散長 ~1 nm → 酸化膜内部へは
   ほぼ入らず、界面第 1 セルのみ平衡化する。テストはこれを前提とする)。
   新しいアクセサを `materials.hpp/cpp` に追加:
   ```cpp
   double segregation_m(const Dopant& d, double temp_k);     // seg_m0*exp(-seg_e/kT)
   double segregation_h(const Dopant& d, double temp_k);     // seg_h0*exp(-seg_he/kT)
   double oxide_diffusivity(const Dopant& d, double temp_k); // dox0*exp(-eox/kT); 0 if dox0<=0
   ```
3. **`DiffusionSolver` コンストラクタの変更**
   (`include/cprocess/diffusion.hpp`):
   ```cpp
   // Per-cell material id: 0 = silicon (full Fair model), 1 = oxide
   // (constant D_ox + segregation exchange with Si), 2 = other (frozen).
   // Empty vector (default) derives ids from solve_mask: mask=1 -> 0,
   // mask=0 -> 2, i.e. the pre-P1-4 behavior.
   DiffusionSolver(const Mesh& mesh, std::vector<char> solve_mask,
                   std::ostream* log = nullptr,
                   std::vector<int> cell_mat = {});
   ```
   メンバ `std::vector<int> mat_;` を追加。コンストラクタで空なら
   mask から導出。`mask_` の意味を「解に参加するセル (mat 0 または 1)」に
   更新する: `mask_[i] = (mat_[i] == 0 || mat_[i] == 1)`。
   (呼び出し側互換のため、`cell_mat` 非空時は `solve_mask` は無視して
   mat から mask を再導出する。)
4. **`build()` の変更** (`src/diffusion.cpp`):
   - `FaceKind` に `kSegregation = 4` を追加。面判定を材料ベースに変更:
     両側 active かつ同材料 → `kInternal`; 片側 Si (mat 0)・片側 oxide
     (mat 1) → `kSegregation`; active/frozen 境界・外部境界 → 従来どおり
     `kBoundOwner`/`kBoundNeigh`/`kInactive`。
   - `FGeom` に `double area = 0;` を追加し、kSegregation 面で
     `area = norm(f.S)` を保存。
   - CSR パターン構築の隣接対追加ループと `fslot_` 構築を
     `kind == kInternal || kind == kSegregation` に拡張する
     (界面の交差項 (iS,iO), (iO,iS) のスロットが必要)。
   - 面彩色: `writes_cells` で kSegregation を kInternal と同様に
     owner+neigh 両セルへ書く面として扱い、active_faces に含める。
   - `bool has_segregation_ = false;` を kSegregation 面が 1 つでも
     あれば true にする。
5. **`assemble()` の変更**: シグネチャに界面係数を渡す。最小変更として
   引数 `double h_seg, double m_seg` を追加する
   (`assemble(dcell, cold, bcface, cgrad, dt, reaction, nonortho,
   h_seg, m_seg, rhs, grad)`; 偏析なしの呼び出しは `h_seg = 0`)。
   面ループに分岐を追加。Si セルを `iS`、oxide セルを `iO`
   (owner/neigh のどちらが Si かは `mat_` で判定)、`A_f = g.area` として
   後退オイラー方程式に界面項を陰的に入れる:
   ```
   Si 行:  (V_S/dt)·C_S + ... + h·A_f·C_S − h·A_f·m·C_O = (V_S/dt)·C_S_old
   Ox 行:  (V_O/dt)·C_O + ... + h·A_f·m·C_O − h·A_f·C_S = (V_O/dt)·C_O_old
   ```
   実装 (交差スロットは fslot_ から):
   ```cpp
   } else if (g.kind == kSegregation && h_seg > 0) {
     const double hA = h_seg * g.area;
     // owner==Si か neigh==Si かで対称に扱う
     A_.val[diag_[iS]]       += hA;
     A_.val[slot(iS, iO)]    -= hA * m_seg;
     A_.val[diag_[iO]]       += hA * m_seg;
     A_.val[slot(iO, iS)]    -= hA;
   }
   ```
   `slot(iS,iO)` は `fslot_[fi][0/1]` を owner/neigh の向きに合わせて選ぶ。
   この項の行和は列ごとに符号が打ち消すので**総量 (Σ C·V) は厳密に保存**
   する (dC_S·V_S/dt = +F·A_f, dC_O·V_O/dt = −F·A_f)。
   注意: m ≠ 1 のとき行列は**非対称**になる。
6. **ソルバー切替** (`run()` / `run_ted()` の線形ソルブ箇所):
   `has_segregation_ && h_seg > 0` のとき CG をスキップして最初から
   `bicgstab_ilu0` を使う (CG は対称前提)。それ以外は現状どおり
   CG → BiCGSTAB フォールバック。
7. **`run()` の変更**:
   - 種ループの per-cell 拡散係数: `mat_[i]==0` は現状の Fair +
     field-enh、`mat_[i]==1` は `oxide_diffusivity(dp, o.temp)`
     (field-enh なし)、frozen は 0。
   - 電荷中性 (`nni`) の集計は **Si セルのみ**。oxide セルは `nni=1`。
   - Picard ループ内の `assemble` 呼び出しに
     `h = segregation_h(dp, o.temp)`, `m = segregation_m(dp, o.temp)` を渡す
     (kSegregation 面が無ければ h=0 で従来と同一)。
   - `run_ted()` は同じ assemble シグネチャ変更に追従する。ψ (格子間原子)
     の式は偏析なし (`h_seg = 0`) のまま。
8. **`proc::diffuse` / `proc::diffuse_ted`** (`src/process.cpp`):
   `st.region_material` から cell_mat を構築するヘルパを追加:
   ```cpp
   // 0 = silicon (or untagged), 1 = oxide/sio2, 2 = everything else
   std::vector<int> material_ids(const SimState& st);
   ```
   solver 構築を
   `DiffusionSolver solver(st.mesh, silicon_mask(st), log, material_ids(st));`
   に変更。`silicon_mask` 自体 (implant が使用) は変更しない。
9. proc:: の**新規関数は無い** (diffuse のシグネチャ不変) ため、CLAUDE.md の
   規約上 pybind / Simulation メソッドの追加は不要。Python 側はテストのみ
   追加する (下記)。

## テスト仕様

`tests/test_segregation.cpp` (新規、`CMakeLists.txt` foreach に
`segregation` を追加)。共通セットアップ: 箱メッシュ 0.05×0.05×0.4 µm、
`make_box_mesh(0, 0.05e-4, 0, 0.05e-4, 0, 0.4e-4, 2, 2, 40)`
(z セル高 0.01 µm)。上部 0.1 µm (centroid z > 0.3e-4) を oxide、
残りを Si とする cell_mat を直接組み、`DiffusionSolver` を素で呼ぶ
(proc 経由テストは 4 で実施)。初期値: Si セルに B = 1e18、oxide は 0。

1. **平衡分配比**: B、1100 °C (1373.15 K)、60 min (`dt=10 s`) 拡散後、
   界面直上の oxide 第 1 セル層と界面直下の Si 第 1 セル層の平均濃度比
   `C_ox / C_si` が `1 / m(1373 K)` = 1/0.369 ≈ **2.71** と相対差 < 10%。
   (D_ox が極小のため oxide 第 1 セルのみ平衡化する — 比較は必ず界面
   隣接セルで行う。)
2. **総量保存**: 同ケースで Σ C·V (Si + oxide 全セル) の前後相対差
   **< 0.5%** (界面項は保存形なので実際は線形ソルバー許容誤差程度)。
3. **B の dose loss**: 1000 °C (1273.15 K)、30 min。Si セルの Σ C·V が
   初期値から **5% を超えて減少**する (期待値 ≈ 10%: 界面 oxide 1 セル
   (0.01 µm) が C_si/0.30 まで満ちるため ≈ (0.01/0.30)/0.30 ≈ 11%)。
4. **P は Si に残る**: 同条件で B の代わりに P。Si の Σ C·V の減少が
   **< 1%** (期待値 ≈ 0.33%: 界面 oxide 1 セル × 1/m = 1/10)。
5. **回帰 (Si 単材料)**: cell_mat を全 0 にした solver と、従来
   コンストラクタ (mask 全 1) の solver で同一の B 拡散
   (1000 °C, 10 min) を実行し、全セル濃度の最大相対差 < 1e-12
   (kSegregation 面が無ければコードパスは同一)。
6. **非対称ソルバー経路**: テスト 1 の実行で例外が出ないこと
   (bicgstab_ilu0 経路の確認; 収束失敗の throw が無ければ PASS)。

Python (`python/test_comprehensive.py` に追加、P1-6 完了後):

1. `Simulation().mesh(0.2, 0.2, 0.4, 4, 4, 40).init("B", 1e18)
   .oxidize(30, 1000).diffuse(30, 1000)` を実行。
   `sim.field("B")` と `sim.cell_centroids` から、酸化前の Si 領域の
   B 総量 (Σ C·V; 体積は等分割セルとして計算してよい) がフロー後に
   5% 以上減少していること、かつ全 B 総量 (Si+oxide) は diffuse 前後で
   1% 以内で保存されていること。

既存テスト: `ctest --test-dir build` の全テスト (特に diffusion / ted /
integration) が PASS すること。assemble のシグネチャ変更は内部 private
なので外部影響なし。

## 完了条件 (DoD)

- [ ] `Dopant` に seg_m0/seg_e/seg_h0/seg_he/dox0/eox + 上表の既定値、
      `segregation_m/h`, `oxide_diffusivity` アクセサ追加
- [ ] `DiffusionSolver` コンストラクタに `cell_mat` (省略時は従来挙動と
      bit 互換) 追加、`kSegregation` 面種 + CSR/彩色/fslot 拡張
- [ ] 界面項が陰的 (両行の対角 + 交差項) に入り、m≠1 で bicgstab_ilu0 に
      切り替わる
- [ ] `proc::diffuse`/`diffuse_ted` が `material_ids(st)` を渡す
- [ ] `tests/test_segregation.cpp` の 6 テスト PASS、foreach 登録
- [ ] Python テスト (oxidize + diffuse で B dose loss) PASS
- [ ] 既存全テスト PASS、コミットメッセージに `P1-4` を含める

## やらないこと

- SiO₂ 以外の材料 (nitride/poly) の拡散・偏析 — 材料テーブルの一般化は
  P1-9。本タスクの cell_mat は 0/1/2 の 3 値のみ
- 酸化**中**の移動界面での snow-plow / dose loss の時間連成 (P2-4)。
  本タスクは静止界面の偏析のみ
- 界面トラップ (3 相目の界面固有サイト) の明示的な面密度状態 — 2 相
  交換 + 平衡係数で近似する
- 偏析係数の濃度依存・応力依存 (P3)
- パラメータのユーザ上書き機構 (P1-10)
