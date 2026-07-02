# P1-9: 多材料拡散 (SiO₂ / Si₃N₄ / poly 中の拡散)

**依存: P1-4 (界面偏析)。** `DiffusionSolver` の cell_mat 引数、
`kSegregation` 面種、非対称ソルバー切替、`Dopant::dox0/eox` と
`segregation_m/h` の機構をそのまま一般化する。

## 目的

P1-4 完了時点で解に参加する非 Si 材料は oxide のみ (cell_mat 0/1/2)。
本タスクでバイナリの solve-mask 的な扱いを**セルごとの材料 id** に置き換え、
nitride (実質バリア)・poly (粒界拡散による高速経路)・gas (凍結) を含む
材料依存拡散を実装する。Fair モデル / 電場増速 / TED 増速は Si セル限定、
他材料は単純 Fick 拡散。偏析界面 (P1-4 の機構) を全材料ペア境界に適用する
(既定の偏析係数は 1、Si/SiO₂ ペアのみ P1-4 の値)。

## 現状コード

- `src/diffusion.cpp` / `include/cprocess/diffusion.hpp` — P1-4 後:
  `DiffusionSolver(mesh, solve_mask, log, cell_mat)` (cell_mat: 0=Si,
  1=oxide, 2=frozen)、`kSegregation` 面、`assemble(..., h_seg, m_seg, ...)`、
  `has_segregation_` → bicgstab_ilu0 切替、`material_ids(st)` ヘルパ
  (process.cpp)。
- `src/materials.cpp` / `include/cprocess/materials.hpp` — `Dopant` に
  Fair 係数 + P1-4 の `dox0/eox`, `seg_*`。`oxide_diffusivity` アクセサ。
- `src/process.cpp` — `silicon_mask(st)` (implant 用、変更しない)、
  `set_region` の許可材料名: silicon/si/oxide/nitride/poly/polysilicon/gas。
  `deposit` の許可材料名に加えて sio2/si3n4 の別名あり。

## 実装手順

1. **材料 id enum** (`include/cprocess/materials.hpp`):
   ```cpp
   enum MatId : int { kMatSi = 0, kMatOxide = 1, kMatNitride = 2,
                      kMatPoly = 3, kMatGas = 4 };
   // Name -> MatId. "silicon"/"si" -> kMatSi, "oxide"/"sio2" -> kMatOxide,
   // "nitride"/"si3n4" -> kMatNitride, "poly"/"polysilicon" -> kMatPoly,
   // "gas" and anything unknown -> kMatGas.
   MatId material_id(const std::string& name);
   ```
   P1-4 の cell_mat 値 0/1/2 とは、0/1 が同一・2 (frozen) が kMatGas=4 に
   相当する形で互換を取る (下記 3)。
2. **`Dopant` の材料別拡散係数** (`include/cprocess/materials.hpp`):
   ```cpp
   // ── Diffusivity outside silicon (P1-9): plain Fickian D0*exp(-E/kT) ──
   // (dox0/eox は P1-4 で追加済み — SiO2 用としてそのまま使う)
   double dnit0 = 0.0, enit = 0.0;    // Si3N4: 0 => perfect barrier
   double dpoly0 = 0.0, epoly = 0.0;  // poly-Si: GB-enhanced, ~10x cryst. Si
   ```
   `kDopants` 既定値 (教科書レベル、較正可。poly は 1000 °C で結晶 Si の
   intrinsic 拡散係数の約 10 倍になるよう活性化エネルギーを Si の主項に
   合わせて前指数を 10 倍相当にしたもの):
   | 種 | dnit0 | enit | dpoly0 (cm²/s) | epoly (eV) | 参考 D_poly(1273 K) |
   |----|-------|------|----------------|-----------|----------------------|
   | B  | 0     | 0    | 7.6  | 3.46 | 1.5e-13 (Si intrinsic ≈ 1.5e-14) |
   | P  | 0     | 0    | 40.0 | 3.66 | 1.3e-13 (Si intrinsic ≈ 1.3e-14) |
   | As | 0     | 0    | 1.1  | 3.44 | 2.6e-14 (Si intrinsic ≈ 2.7e-15) |
   | Sb | 0     | 0    | 5.3  | 3.65 | 1.9e-14 (Si intrinsic ≈ 1.8e-15) |
   汎用アクセサを追加し、P1-4 の `oxide_diffusivity` はこれに委譲:
   ```cpp
   // D in material `mat` at temp_k for n/ni = nni.
   // kMatSi -> full Fair model dopant_diffusivity(d, temp_k, nni);
   // kMatOxide/kMatNitride/kMatPoly -> plain Arrhenius from the pairs above
   // (nni ignored); kMatGas -> 0.
   double material_diffusivity(const Dopant& d, MatId mat, double temp_k,
                               double nni);
   ```
3. **`DiffusionSolver` の一般化** (`src/diffusion.cpp`):
   - `mat_` の値域を MatId (0..4) に拡張。「解に参加」の定義:
     `active(i) = (material_diffusivity(…) が正になり得る材料)` ではなく
     単純に `mat_[i] != kMatGas`。ただし、その種の D が 0 の材料
     (既定の nitride) のセルは**その種のソルブでは凍結扱い**にする
     (下記の面規則参照)。実装: 種ごとの dcell 計算で D=0 のセルは
     dcell=0 とし、面種は build() で材料ペアから静的に決めるが、
     assemble 時に dcell が両側 >0 の面のみフラックスを立てる
     (kInternal は既存の `dP<=0 || dN<=0 → continue` ガードが既にこれを
     行う; kSegregation にも同じガードを追加する)。
   - 面種の規則 (build()):
     - 両側 gas → `kInactive`; 片側 gas → gas でない側の
       `kBoundOwner`/`kBoundNeigh` (従来の frozen 境界と同じ)
     - 両側同材料 (非 gas) → `kInternal`
     - 異材料ペア (非 gas 同士) → `kSegregation`
   - 後方互換: `cell_mat` が空のときは従来どおり solve_mask から導出
     (mask=1 → kMatSi, mask=0 → kMatGas)。P1-4 の 3 値 (0/1/2) 呼び出しは
     2 → kMatGas と読み替えても意味が変わらない (2 と 4 を同義とする)。
   - `run()`/`run_ted()` の per-cell 係数:
     `dcell[s][i] = material_diffusivity(dp, mat_[i], o.temp, nni[i])`;
     電場増速 (`o.field_enh`) の分岐と TED のペア増速
     `(1-fi)+fi*S[i]` は **`mat_[i] == kMatSi` のときだけ**掛ける。
     電荷中性 `nni` の集計も Si セルのみ (P1-4 と同じ)、非 Si は nni=1。
     TED の格子間原子 ψ の拡散係数 `dI[i]` も Si セルのみ非零 (現状の
     mask ベースから mat ベースへ書き換え)。
   - **偏析係数の材料ペア一般化**: 界面係数 h, m は「Si↔SiO₂ ペア」のみ
     `segregation_m/h(dp, T)` を使い、**それ以外の全ペアは m = 1、
     h = segregation_h(dp, T)** (輸送は速く、分配は等分) とする。
     実装: assemble に h_seg/m_seg のスカラーを渡す現行 (P1-4) 形を、
     面ごとに引けるよう `std::array<double, 2>` を面種別に持つのではなく、
     **assemble の引数を「材料ペア (min(matP,matN), max(matP,matN)) →
     (h, m) の 5×5 小テーブル `SegTable`」に変更**する:
     ```cpp
     struct SegTable { double h[5][5] = {}; double m[5][5]; };  // m init to 1
     ```
     run() が種・温度ごとに SegTable を埋めて assemble に渡す。
     m の向き規約は P1-4 と同じ (`m = C_lo/C_hi` ではなく
     **m = C_{mat 小さい側 id} / C_{mat 大きい側 id}**、つまり
     Si(0)/oxide(1) ペアでは m = C_si/C_ox — P1-4 の実装と一致)。
     m=1 のペアのみなら行列は対称のままなので、`bicgstab_ilu0` への切替
     条件を「m ≠ 1 の kSegregation 面が存在する」に精密化する
     (Si/SiO₂ 界面が無ければ CG のまま)。
4. **`proc::diffuse` / `proc::diffuse_ted`** (`src/process.cpp`):
   P1-4 の `material_ids(st)` を `material_id(材料名)` ベースに書き換え、
   untagged セルは kMatSi。`silicon_mask` は implant 用に**そのまま残す**。
5. proc:: の**新規関数は無い** (diffuse シグネチャ不変) ため pybind /
   Simulation メソッドの追加は不要 (CLAUDE.md 規約)。Python はテストのみ。

## 実装順の注意

build() の面種規則・SegTable 化は P1-4 のコードを置き換える差分になる。
P1-4 のテスト (test_segregation.cpp) は一切変更せずに通ること。

## テスト仕様

`tests/test_multimaterial.cpp` (新規、`CMakeLists.txt` foreach に
`multimaterial` を追加)。

1. **nitride バリア**: 箱メッシュ 0.05×0.05×0.4 µm (2×2×40、セル高
   0.01 µm)、上部 0.1 µm を kMatNitride、下部を kMatSi とする cell_mat で
   solver を直接構築。Si 側に B の Gaussian (Rp=0.2 µm 相当の深さ、
   ピーク 1e19) を置き 1000 °C 30 min 拡散。拡散後の nitride セル内の
   B 総量が全 B 総量の **< 0.1%** (D_nit=0 かつ D=0 側とのフラックス
   ガードにより厳密に 0 のはず)。かつ全 B 総量の保存誤差 < 0.5%。
2. **poly 高速経路**: 同一形状の 1D カラム (2×2×40、0.4 µm) を 2 本、
   一方は全セル kMatSi、他方は全セル kMatPoly。同一の初期 Gaussian B
   (中心 0.2 µm、σ=0.02 µm、ピーク 1e18 — intrinsic 条件で電場/濃度依存の
   影響を避ける)、1000 °C 30 min (`field_enh=false` で実行)。
   深さプロファイルの分散 σ² (セル centroid z の濃度重み付き分散) の増分比
   `(σ²_poly − σ²_0) / (σ²_si − σ²_0)` が **[5, 20] の範囲**
   (期待値 = D_poly/D_si ≈ 10)。
3. **Si 単材料の回帰**: 全セル kMatSi の cell_mat 付き solver と、
   cell_mat 省略 (mask 全 1) の solver で同一の B+P 連成拡散
   (1000 °C 10 min、field_enh=true) を実行し、各種の総 dose とピーク値の
   相対差 **< 1e-6** (コードパスがほぼ同一なので実際は ~0; bit 一致は
   要求しない)。
4. **異材料ペアの総量保存**: Si(下 0.2 µm) / poly(中 0.1 µm) /
   oxide(上 0.1 µm) の 3 層カラムに B 一様 1e18 (Si と poly のみ)、
   1000 °C 30 min。全 B 総量の保存誤差 < 0.5%。Si/poly 界面 (m=1) を
   はさんだ隣接セル濃度比が 1 ± 10% (等分配)。

既存テスト: `ctest` の diffusion / ted / segregation / integration を含む
全テストが**変更なしで** PASS すること (Si-only メッシュの挙動は
テスト 3 の基準で同一性を担保)。

Python (`python/test_comprehensive.py` に追加):

1. `Simulation().mesh(0.2, 0.2, 0.4, 4, 4, 40)` →
   `implant` (B, Gaussian, rp 指定で浅め) → `deposit("nitride", 0.1)` →
   `diffuse(30, 1000)`。`sim.field("B")` と `sim.cell_centroids` から、
   nitride 層 (z > 0.4 µm、deposit で上に載る) 内の B 総量が全体の
   < 0.1% であること。
2. 既存 Python テストが全て PASS。

## 完了条件 (DoD)

- [ ] `MatId` enum + `material_id()` + `material_diffusivity()` 追加、
      `Dopant` に dnit0/enit/dpoly0/epoly + 上表の既定値
- [ ] `DiffusionSolver` が MatId ベース (gas 凍結、材料ペア → kSegregation、
      D=0 側があるペアはフラックスなし)、SegTable による面ごとの h/m
- [ ] Fair / 電場増速 / TED 増速 / nni 集計 / dI が Si セル限定
- [ ] m=1 のみなら CG のまま、m≠1 ペア存在時のみ BiCGSTAB
- [ ] `material_ids(st)` が 5 材料を返す。`silicon_mask` は不変
- [ ] `tests/test_multimaterial.cpp` の 4 テスト PASS、foreach 登録、
      test_segregation.cpp は無変更で PASS
- [ ] Python テスト (deposit("nitride") + diffuse) PASS、既存全テスト PASS
- [ ] コミットメッセージに `P1-9` を含める

## やらないこと

- 材料別の完全な較正 DB (Advanced Calibration 相当) — 上表は教科書レベルの
  既定値。ユーザ上書きは P1-10
- resist 中の拡散 (resist は stack 側の概念で拡散メッシュには現れない)
- In/C/F/Ge 等の新ドーパント (P2-8)、SiGe (P3)
- poly の粒界/粒内 2 経路モデル — 単一の実効 Arrhenius で近似
- 非 Si 材料での活性化・クラスタリング (P1-3/P2-2)
- 界面での材料ペア依存 h の較正 (既定は全ペア共通の h(T))
