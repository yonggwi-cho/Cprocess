# Sentaurus Process リファレンス（機能・モデル・数値計算手法）

> Synopsys **Sentaurus Process (sprocess)** の機能一覧・物理モデル一覧・
> 数値計算手法を、可能な限り数式付きで整理したリファレンス。
> Cprocess を業界標準へ寄せるための実装目標として用いる。
>
> 作成日: 2026-06-15 / 出典は末尾「参照」を参照。
> 注: 数式は公開チュートリアル・学位論文・標準教科書に基づく標準形であり、
> Synopsys の内部実装係数そのものではない（校正値は `pdbSet` で設定される）。

---

## 1. 機能一覧

| カテゴリ | 機能 | コマンド例 |
|---|---|---|
| 構造定義 | 1D/2D/3D メッシュ、領域・材料定義 | `line`, `region`, `init` |
| メッシュ | 適応細分化、リメッシュ、変換 | `refinebox`, `grid remesh`, `transform` |
| イオン注入 | 解析 (Pearson) / MC (Crystal-TRIM) | `implant` |
| 拡散・アニール | 平衡〜非平衡（five-stream）、酸化雰囲気 | `diffuse` |
| 熱酸化 | Deal-Grove + 応力 + 移動境界 | `diffuse ... O2/H2O` |
| デポジション | 等方/異方、材料積層 | `deposit` |
| エッチング | 等方/異方、マスク | `etch`, `mask`, `photo`, `strip` |
| 活性化 | クラスタ/析出/固溶 | `pdbSet ... ActiveModel` |
| 応力 | 粘弾性/塑性、格子不整合 | `stress`, `pdbSet` |
| パラメータDB | 階層 KV の実行時上書き | `pdbSet`, `pdbGet`, `pdbSetDouble` |
| 出力 | TDR 構造、1D プロファイル抽出 | `struct tdr=`, `WritePlx`, `SetPlxList` |
| SiC | 2H/4H/6H/3C ポリタイプ、ミスカット | `set4H-SiC`, `init wafer.orient=` |

### 座標系規約
- 1D: x はウェハ内部（下向き）が正。2D もこれを踏襲し、y が横方向。

---

## 2. 物理モデル一覧

### 2.1 イオン注入

#### (a) 解析モデル — Pearson-IV / dual-Pearson
深さ分布 f(x) を 4 モーメント（飛程 Rp, 標準偏差 σ=ΔRp, 歪度 γ, 尖度 β）で表す。
Pearson 系は次の微分方程式の解として定義される:

```
df/dx = (x − a) f / (b0 + b1 (x − a) + b2 (x − a)²)
```

係数はモーメントから:

```
a  = −σ γ (β + 3) / D
b0 = −σ² (4β − 3γ²) / D
b1 = a
b2 = −(2β − 3γ² − 6) / D
D  = 10β − 12γ² − 18
```

チャネリングテールは 2 つの Pearson を重み付き加算（dual-Pearson）:
`f(x) = (1−r) f1(x) + r f2(x)`。

横方向は通常ガウシアン:
`F(x,y) = f(x) · (1/(√(2π) σ_lat)) exp(−y²/(2 σ_lat²))`。

ドーズ整合: `∫ f(x) dx = Φ`（ドーズ Φ [cm⁻²]）。
校正は **Sentaurus Calibration Library**（Varian 等の実測モーメント表）。

#### (b) Monte Carlo — Crystal-TRIM (BCA)
二体衝突近似（Binary Collision Approximation）。散乱角は核衝突の散乱積分:

```
θ: cos(θ/2) から、衝突径数 p と ZBL 万有遮蔽ポテンシャル
V(r) = (Z1 Z2 e²)/r · Φ(r/a),  a = 0.8854 a0 / (Z1^0.23 + Z2^0.23)
Φ(x) = Σ ci exp(−di x)   (ZBL: ci={0.1818,0.5099,0.2802,0.02817},
                                  di={3.2,0.9423,0.4029,0.2016})
```

電子的阻止能（Lindhard-Scharff、低速領域）:
`Se(E) = k_LS √E`。
結晶ターゲットではチャネリング（連続ストリングポテンシャル）。
SiC では `LSS.pre`, `nloc.pre`, `nloc.exp`, `surv.rat` で校正。

#### (c) 注入ダメージ
Kinchin-Pease の変位数:

```
N_d(T) = 0  (T < Ed),  = 1  (Ed ≤ T < 2Ed/0.8),
       = 0.8 T / (2 Ed)  (T ≥ 2Ed/0.8)
```

(Ed: 変位閾エネルギー、Si で約 15 eV)。"+1" モデルで残留格子間原子を後段 TED の
初期条件に与える。蓄積ダメージで動的非晶質化。

### 2.2 拡散モデル階層
`pdbSet Silicon <Species> DiffModel <model>` で選択。

#### (a) Fermi（平衡・点欠陥は平衡）
実効拡散係数を Fermi 準位（n/ni）の多項式で:

```
D_A(n) = D⁰ + D⁻ (n/ni) + D⁼ (n/ni)² + D⁺ (ni/n)
各成分: D^c = D0^c · exp(−E^c / kT)   (Arrhenius)
```

拡散方程式（電界ドリフト込み）:

```
∂C_A/∂t = ∇·[ D_A (∇C_A + Z C_A (∇n/n)) ]     (Z: 電荷符号)
```

内部電界増速因子 h:
`h = 1 + (Nnet/2ni) / √((Nnet/2ni)² + 1)`、`n/ni = (Nnet/2ni) + √(...)`。

#### (b) ChargedFermi
各荷電状態を陽に区別した Fermi モデル（上記の荷電別総和）。

#### (c) Pair / ChargedPair
ドーパント-欠陥ペア (AI), (AV) を介した拡散。点欠陥は準平衡。
有効拡散係数は点欠陥過飽和 (C_I/C_I*) に比例:

```
D_A^eff = D_A^I (C_I / C_I*) + D_A^V (C_V / C_V*)
```

#### (d) React（five-stream / Dunham モデル）★最重要
5 連立: ドーパント A⁺、格子間 I、空孔 V、ペア AI、ペア AV。

**反応:**
```
A⁺ + Iⁱ ⟺ (AI)^(i+1)
A⁺ + Vⁱ ⟺ (AV)^(i+1)
Iⁱ + Vʲ ⟺ 0           (Frenkel 再結合)
(AI)ⁱ + Vʲ ⟺ A⁺
(AV)ⁱ + Iʲ ⟺ A⁺
(AI)ⁱ + (AV)ʲ ⟺ 2A⁺
```

**連続の式:**
```
∂C_A /∂t = −R_AI − R_AV + R_(AI+V) + R_(AV+I) + 2 R_(AI+AV)
∂C_I /∂t = −∇·J_I  + R_AI − R_(I+V) − R_(AV+I)
∂C_V /∂t = −∇·J_V  + R_AV − R_(I+V) − R_(AI+V)
∂C_AI/∂t = −∇·J_AI + R_AI − R_(AI+V) − R_(AI+AV)
∂C_AV/∂t = −∇·J_AV + R_AV − R_(AV+I) − R_(AI+AV)
```

**フラックス（点欠陥・ペア）:**
```
J_I  = −[ Σᵢ D_Iⁱ K_Iⁱ (n/ni)ⁱ ] ∇C_I⁰
J_V  = −[ Σᵢ D_Vⁱ K_Vⁱ (n/ni)ⁱ ] ∇C_V⁰
J_AI = −[ Σᵢ D_(AI)^(i+1) K_AIⁱ (n/ni)ⁱ ] (∇C_AI⁺ + C_AI⁺ (ni/n) ∇(n/ni))
J_AV = −[ Σᵢ D_(AV)^(i+1) K_AVⁱ (n/ni)ⁱ ] (∇C_AV⁺ + C_AV⁺ (ni/n) ∇(n/ni))
```

**反応速度（質量作用則、例）:**
```
R_AI = [Σᵢ k_AIⁱ K_Iⁱ (n/ni)ⁱ] · (C_A⁺ C_I⁰ − C_AI⁺ / K_(A+I)⁰)
R_AV = [Σᵢ k_AVⁱ K_Vⁱ (n/ni)ⁱ] · (C_A⁺ C_V⁰ − C_AV⁺ / K_(A+V)⁰)
R_(I+V) = [Σ_ij k^ij K_Iⁱ K_Vʲ (n/ni)^(i+j)] · (C_I⁰ C_V⁰ − C_I⁰* C_V⁰*)
```

（C_X*: 平衡濃度、K_X: 荷電状態平衡定数、k: 前方反応速度）。
この非平衡結合系が **TED（過渡増速拡散）** を再現する本体。

### 2.3 欠陥・クラスタ・活性化
- {311} 欠陥、SMIC（小格子間クラスタ）、dopant-defect クラスタ。
- `ActiveModel`: `None / Solid / Transient / Precipitation / Cluster / BIC /
  ChargedCluster / ComplexCluster / Equilibrium`。
- 固溶度: `C_ss(T) = C0 exp(−E_ss / kT)`。可動濃度はこれでクランプ、超過分は
  クラスタ/析出相へ。

### 2.4 点欠陥の界面生成・消滅
界面再結合（例 Oxide_Silicon 界面の空孔）:

```
K_surf = A exp(−E / kT)    [cm/s]
J·n = K_surf (C − C*)      (界面フラックス境界条件)
```

### 2.5 熱酸化 — Deal-Grove
酸化膜厚 x_ox の成長:

```
x_ox² + A x_ox = B (t + τ)
線形速度定数: B/A = (B/A)0 exp(−E1/kT)
放物速度定数: B   = B0   exp(−E2/kT)
解: x_ox = (A/2)[ √(1 + (t+τ)/(A²/4B)) − 1 ]
```

雰囲気 O₂（乾式）/H₂O（湿式）で B,A が異なる。圧力依存、結晶方位依存、
酸化増速拡散（OED: 過剰 I 注入）、移動境界（体積膨張比 ~2.2）、粘弾性応力と連成。

### 2.6 応力
粘弾性（Maxwell）構成則 + 格子不整合・熱膨張差・酸化体積膨張を駆動力とする。
応力は拡散係数・偏析・酸化速度に影響。

### 2.7 境界条件
- **HomNeumann**（既定）: 無フラックス `∂C/∂n = 0`。
- **Dirichlet**: 界面平衡濃度を固定（1e10〜1e19 cm⁻³）。
- **偏析**: 界面で `C_1/C_2 = m`（偏析係数、Arrhenius `[Arr 0.75 −0.336]` 等）。

---

## 3. 数値計算手法

| 項目 | Sentaurus の手法 |
|---|---|
| 空間離散化 | ボックス法（有限体積）/ 有限要素、非構造メッシュ |
| メッシュ適応 | 誤差指標ベースの細分化・粗化、移動境界対応 |
| 時間積分 | 後退オイラー / TR-BDF2（陰、A-安定） |
| 非線形 | **Newton-Raphson**（全種・全欠陥を完全結合、解析ヤコビアン） |
| 線形ソルバ | 疎直接法（PARDISO 系 LU）/ 反復法（ILU 前処理）|
| 並列 | 共有メモリ（`math numThreads=`）、MC は `numThreadsMC=` |
| パラメータ評価 | `[Arrhenius pre Ea]` を実行時に温度で評価 |

### 3.1 離散化（ボックス法）
各制御体積 Ωi で連続の式を積分し発散定理を適用:

```
|Ωi| (C_i^{n+1} − C_i^n)/Δt = Σ_{j∈neigh} A_ij J_ij + |Ωi| R_i
J_ij = D_ij (C_j − C_i)/d_ij + (ドリフト項)
```

### 3.2 Newton 法
非線形残差 F(u)=0（u は全種濃度ベクトル）を:

```
J(u^k) Δu = −F(u^k),   u^{k+1} = u^k + Δu,   J = ∂F/∂u
```

five-stream は強結合のため Newton（Picard 不動点ではなく）が必須。

### 3.3 Arrhenius パラメータ表現
`[Arrhenius A Ea] → A·exp(−Ea/kT)`。`pdbSet`/`pdbGet`/`pdbDelayDouble` で
取得・上書き。例: `pdbSet Silicon Boron Dstar {[expr 4.0*[pdbDelayDouble Silicon Boron Dstar]]}`。

---

## 4. Cprocess とのギャップ（要点）

詳細は `docs/sentaurus_gap_analysis.md` を参照。本書との対応:

- §2.1(a) Pearson-IV → Cprocess はガウシアンのみ（モーメント未対応）。
- §2.2(a) Fermi → Cprocess は **実装済み相当**（中性/+/−/++ の荷電多項式）。
- §2.2(d) five-stream → Cprocess **未実装**（最重要ターゲット）。
- §2.5 Deal-Grove 酸化 → **未実装**。
- §3.2 Newton → Cprocess は Picard 反復（five-stream 移行時に Newton 化が必要）。
- §1 `pdbSet` → Cprocess はハードコード（`materials.cpp`）。

---

## 参照

- Sentaurus Process チュートリアル（コマンド/2D プロセス）
  <https://ghzphy.github.io/Sentaurus_Training/sp/sp_3.html>
- Sentaurus Process チュートリアル（モデル指定 `pdbSet`/Arrhenius）
  <https://ghzphy.github.io/Sentaurus_Training/sp/sp_4.html>
- Sentaurus Process チュートリアル（SiC プロセス）
  <https://ghzphy.github.io/Sentaurus_Training/sp/sp_15.html>
- Five-Stream Dunham Diffusion Model（TU Wien, Hollauer 学位論文）
  <https://www.iue.tuwien.ac.at/phd/hollauer/node24.html>
- Synopsys TCAD / Process Simulation 製品ページ
  <https://www.synopsys.com/manufacturing/tcad/process-simulation.html>
- Sentaurus Process User Guide（Y-2006.06、公開版マニュアル）
  <https://dunham.ece.uw.edu/protected/PDFManual/data/sprocess_ug.pdf>
- B.J. Mulvaney et al. / R.B. Fair（拡散係数）, B.E. Deal & A.S. Grove
  J. Appl. Phys. 36 (1965) 3770（熱酸化）
