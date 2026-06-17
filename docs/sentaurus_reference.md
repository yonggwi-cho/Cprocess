# Sentaurus Process 技術リファレンス完全版
## 機能一覧・モデル一覧・数値計算手法（数式付き）

> 出典: Synopsys TCAD 公開チュートリアル (ghzphy.github.io/Sentaurus_Training/)、
> TU Wien 博士論文 (iue.tuwien.ac.at)、標準半導体プロセス教科書、IEEE 論文。
> Synopsys 社内の未公開係数は含まない。
> 作成日: 2026-06-15

---

## 目次

1. [機能一覧（コマンドリファレンス）](#1-機能一覧)
2. [物理モデル一覧（数式付き）](#2-物理モデル)
   - 2.1 イオン注入
   - 2.2 拡散モデル階層
   - 2.3 熱酸化（Deal-Grove）
   - 2.4 欠陥・ダメージ・TED
   - 2.5 ドーパント活性化・クラスタ
   - 2.6 応力
   - 2.7 界面・境界条件
3. [数値計算手法（数式付き）](#3-数値計算手法)
   - 3.1 空間離散化（ボックス法 FVM）
   - 3.2 時間積分
   - 3.3 非線形ソルバー（Newton-Raphson）
   - 3.4 線形ソルバー
   - 3.5 Arrhenius パラメータ評価
4. [パラメータDB（pdbSet）](#4-パラメータdb)
5. [単位系・座標系](#5-単位系座標系)
6. [参照](#6-参照)

---

## 1. 機能一覧

### 1.1 グリッド・構造定義コマンド

```tcl
# グリッド線定義（1D/2D/3D）
line x|y|z  location= <値><単位>  spacing= <値><単位>  tag= <ラベル>

# 領域・材料定義
region <Material>  xlo= <tag>  xhi= <tag>  [ylo=  yhi=  zlo=  zhi=]

# 基板初期化（ドーピング設定）
init  concentration= <値><cm-3>  field= <元素名>
init  tdr= <ファイル>                          # 構造リストア

# ウェハ仕様（SiC など）
init  wafer.orient= {a1 a2 c}  notch.direction= {a1 a2 c}
      miscut.tilt= <角度>  miscut.toward= {a1 a2 c}
```

### 1.2 メッシュ制御コマンド

```tcl
# リメッシュパラメータ
grid  set.min.normal.size= <値>       # 最小法線方向間隔
      set.normal.growth.ratio.2d= <値>  # メッシュ拡大比

# 適応細分化ボックス
refinebox  <Material>  min= {x y}  max= {x y}  xrefine= {v1 v2 v3}  yrefine= {v1 v2 v3}  add
grid remesh                           # 細分化実行

# 酸化移動境界制御
pdbSet Oxide Grid perp.add.dist <値>  # 酸化界面メッシュ幅（既定 10e-7 cm）
pdbSetBoolean Grid DoNotMove.Reaction 1   # 3D で境界移動なし（高速）
```

### 1.3 プロセスステップコマンド

```tcl
# イオン注入（解析）
implant <Species>  dose= <値><cm-2>  energy= <値><keV>
                   tilt= <度>  rotation= <度>
                   [sentaurus.mc | crystaltrim]      # MC モード切替
                   [particles= <粒子数>]

# 拡散・アニール
diffuse  temperature= <値><C>  time= <値><min>
         [gas_flow= <ID>]          # 酸化雰囲気
         [O2 | H2O | N2]           # 雰囲気短縮形
         [pressure= <atm>]

# 雰囲気ガス定義
gas_flow  name= <ID>  pressure= <atm>  flowO2= <l/min>  flowN2= <l/min>

# デポジション
deposit  material= {<Material>}  type= isotropic|anisotropic
         time= <値>  rate= {<値>}       # または
         thickness= <値><nm>

# エッチング
etch  material= {<Material>}  type= isotropic|anisotropic
      time= <値>  rate= {<値>}
      [mask= <maskname>]  [isotropic.overetch= <値>]

# マスク定義
mask  name= <ID>  left= <値>  right= <値>  [segments= {<値>...}]  [negative]
photo  mask= <ID>  thickness= <値>
strip  <Material>                  # 材料剥離（例: strip PhotoResist）

# 構造保存・リストア
struct  tdr= <stem>                # <stem>_fps.tdr, _bnd.tdr, _msh.cmd を生成
struct  smesh= <stem>
```

### 1.4 後処理・出力コマンド

```tcl
select  z= <式>                    # 測定量選択（例 z=1 で厚さ測定）
layers                             # 層厚・積分値レポート
SetPlxList  {BTotal AsTotal NetActive ...}
WritePlx    <ファイル>  [y= <値>]  [<Material>]
plot.1d     [...]
transform   reflect|rotate|cut     # 構造変換
contact     name= <ID>  ...        # 電極定義（デバイスシミュレーション連携）
print       ...
pdbSet      InfoDefault <0|1|2>    # 出力冗長度
exit | quit
```

### 1.5 SiC 固有コマンド

```tcl
set4H-SiC   # または set2H-SiC, set6H-SiC, set3C-SiC
AdvancedCalibration 4H-SiC
math numThreadsMC= <整数>          # MC 並列スレッド数
math numThreads=  <整数>           # 拡散ソルバー並列数
```

---

## 2. 物理モデル

---

### 2.1 イオン注入

#### (a) 解析モデル — Pearson IV 分布

深さ分布を 4 つの統計モーメントで記述する。分布関数 f(x) は以下の微分方程式の解:

```
df(x)/dx = (x − a) f(x) / (b₀ + b₁(x−a) + b₂(x−a)²)
```

係数はモーメントから決定される:

```
D  = 10β − 12γ² − 18
a  = −σ γ (β + 3) / D
b₀ = −σ² (4β − 3γ²) / D
b₁ = a
b₂ = −(2β − 3γ² − 6) / D
```

統計モーメントの定義:
```
μ₁ = Rp          (平均飛程, projected range)
μ₂ = σ²          (分散, straggling²)
γ  = μ₃ / μ₂^(3/2)   (歪度, skewness)
β  = μ₄ / μ₂²        (尖度, kurtosis)
```

Pearson 型の分類（β, γ の値域で自動選択）:
- **IV 型** (最も一般): 閉形式なし、数値積分で規格化
- **VI 型**: β > 6γ²/5 + 1 の領域
- **正規分布** (Gaussian): γ=0, β=3 の特殊ケース

ドーズ規格化: `∫₋∞^∞ f(x) dx = Φ [cm⁻²]`

横方向はガウシアン:
```
F(x, y) = f(x) · (1 / (√(2π) σ_lat)) · exp(−y² / (2 σ_lat²))
```

#### (b) 双峰 Pearson（Dual-Pearson）— チャネリングテール対応

2 つの Pearson 分布の重み付き加算:

```
C(x) = (1 − r) · f₁(x; Rp₁, σ₁, γ₁, β₁) + r · f₂(x; Rp₂, σ₂, γ₂, β₂)
```

9 パラメータ: P₁=(Rp₁, σ₁, γ₁, β₁), P₂=(Rp₂, σ₂, γ₂, β₂), r（重み）。
- P₁: 表面近傍（非チャネリング成分）
- P₂: 深部チャネリングテール成分

係数は **Sentaurus Calibration Library**（Varian 等実測データ）から供給。

#### (c) MC 注入 — Crystal-TRIM (BCA)

**核散乱**: ZBL 万有遮蔽ポテンシャル

```
V(r) = (Z₁ Z₂ e²) / (4πε₀ r) · Φ_ZBL(r / a)

a = 0.8854 a₀ / (Z₁^0.23 + Z₂^0.23)   (遮蔽長, a₀=Bohr半径)

Φ_ZBL(x) = 0.1818 e^(-3.2x) + 0.5099 e^(-0.9423x)
           + 0.2802 e^(-0.4029x) + 0.02817 e^(-0.2016x)
```

散乱角（重心系）は衝突径数 p からの積分:

```
θ_cm = π − 2p ∫_{r_min}^∞ dr / (r² √(1 − V(r)/E_cm − p²/r²))
```

**電子的阻止能**（Lindhard-Scharff、低速領域）:
```
Se(E) = k_LS · √E,    k_LS = 0.0793 (Z₁^(2/3) Z₂^(1/2)) / (Z₁^(2/3) + Z₂^(2/3))^(3/2)
                              · (A₁ + A₂)^(1/2) / (A₁^(3/2) · A₂^(1/2))   [eV/(10¹⁵ cm⁻²)]
```

高速領域（Bethe-Bloch 近似）:
```
-dE/dx ∝ (Z₁² / E) · Z₂ · ln(4 mₑ E / (mₑ + m₁) · 1/I)
```

**チャネリング（Si 結晶、連続ストリング）**:
```
U(ρ) = Z₁ Z₂ e² / (d√(ρ² + 3u₁²))    (Doyle-Turner 形式)
```
ここで d = チェーン間隔, ρ = 弦からの距離, u₁ = 熱振動振幅。
非晶質化が進むにつれてポテンシャルを低減: `U_eff = U_max · (1 − f_amor)`。

---

### 2.2 拡散モデル階層

`pdbSet Silicon <Species> DiffModel <model>` で選択。

#### (a) Fermi モデル（平衡・点欠陥平衡仮定）

拡散フラックス（電界ドリフト込み）:
```
J = −D_A(n) [∇C_A + Z · C_A · (∇n / n)]
```

実効拡散係数（荷電ペアの加算）:
```
D_A(n) = D⁰ + D⁻(n/nᵢ) + D⁼(n/nᵢ)² + D⁺(nᵢ/n)
```

各成分は Arrhenius:
```
D^c = D₀^c · exp(−Eₐ^c / kT)
```

連続の式:
```
∂C_A/∂t = ∇ · [D_A(n) (∇C_A + Z C_A ∇n/n)]
```

電子濃度 n は電荷中性から（全ドーパント合算）:
```
(n − p) = Σ_donors C_D − Σ_acceptors C_A,    n·p = nᵢ²
n/nᵢ = N_net/(2nᵢ) + √((N_net/(2nᵢ))² + 1)
```

内部電界増速因子 h（電界ドリフトを含む実効倍率）:
```
h = 1 + |N_net / 2nᵢ| / √((N_net / 2nᵢ)² + 1)
```

真性キャリア濃度（Si）:
```
nᵢ(T) = 3.87 × 10¹⁶ · T^1.5 · exp(−0.605 / (kT))   [cm⁻³, T in K]
```

**Fair の真空モデル係数（SUPREM 標準値）:**

| 種 | D⁰ [cm²/s] | E⁰ [eV] | D⁻ | E⁻ | D⁼ | E⁼ | D⁺ | E⁺ |
|---|---|---|---|---|---|---|---|---|
| B | 0.037 | 3.46 | — | — | — | — | 0.72 | 3.46 |
| P | 3.85 | 3.66 | 4.44 | 4.00 | 44.2 | 4.37 | — | — |
| As | 0.066 | 3.44 | 12.0 | 4.05 | — | — | — | — |
| Sb | 0.214 | 3.65 | 15.0 | 4.08 | — | — | — | — |

固溶度（Arrhenius）:
```
C_ss(T) = C₀ · exp(−Eₛₛ / kT)
```

| 種 | C₀ [cm⁻³] | Eₛₛ [eV] |
|---|---|---|
| B | 9.25 × 10²² | 0.73 |
| P | 2.45 × 10²³ | 0.62 |
| As | 1.3 × 10²³ | 0.66 |
| Sb | 3.8 × 10²¹ | 0.56 |

#### (b) ChargedFermi モデル

各荷電状態を明示的に合算（Fermi の拡張版）。実効 D は上記と同形だが、
荷電状態ごとの化学ポテンシャルを厳密に扱う。

#### (c) Pair モデル（準平衡点欠陥）

ドーパント A がインタースティシャル I または空孔 V とペアを組んで拡散:
```
D_A^eff = D_AI · (C_I / C_I*) + D_AV · (C_V / C_V*)
```

- `C_I*, C_V*`: 点欠陥の熱平衡濃度
- `C_I / C_I*`: 注入・酸化由来の過飽和率（点欠陥は別方程式で求める）

#### (d) ChargedPair モデル

Pair の拡張。荷電状態別ペアを陽に区別し Fermi 準位依存性を含む。

#### (e) React モデル（Five-Stream / Dunham モデル）★ 最重要

5 種を連立 PDE で結合し TED を再現する。

**反応（イオン化種と荷電欠陥の全結合）:**
```
A⁺ + Iⁱ  ⇌  (AI)^(i+1)           ドーパント-格子間ペア生成
A⁺ + Vⁱ  ⇌  (AV)^(i+1)           ドーパント-空孔ペア生成
Iⁱ + Vʲ  ⇌  0  + (-i-j)e⁻        Frenkel 再結合
(AI)ⁱ + Vʲ ⇌ A⁺ + (1-i-j)e⁻     ペア-空孔消滅
(AV)ⁱ + Iʲ ⇌ A⁺ + (1-i-j)e⁻     ペア-格子間消滅
(AI)ⁱ +(AV)ʲ⇌ 2A⁺ + (2-i-j)e⁻   ペア間消滅
```

**5 連立連続の式:**
```
∂C_{A⁺}/∂t  =  −∇J_A  − R_{AI} − R_{AV} + R_{AI+V} + R_{AV+I} + 2R_{AI+AV}
∂C_I/∂t     =  −∇J_I            − R_{AI} − R_{I+V}  − R_{AV+I}
∂C_V/∂t     =  −∇J_V            − R_{AV} − R_{I+V}  − R_{AI+V}
∂C_{AI}/∂t  =  −∇J_{AI}         + R_{AI} − R_{AI+V} − R_{AI+AV}
∂C_{AV}/∂t  =  −∇J_{AV}         + R_{AV} − R_{AV+I} − R_{AI+AV}
```

**フラックス式:**
```
J_I  = −[Σᵢ D_{Iⁱ} K_Iⁱ (n/nᵢ)ⁱ] ∇C_{I⁰}

J_V  = −[Σᵢ D_{Vⁱ} K_Vⁱ (n/nᵢ)ⁱ] ∇C_{V⁰}

J_{AI} = −[Σᵢ D_{(AI)^{i+1}} K_{AI}ⁱ (n/nᵢ)ⁱ]
           · (∇C_{(AI)⁺} + C_{(AI)⁺} (nᵢ/n) ∇(n/nᵢ))

J_{AV} = −[Σᵢ D_{(AV)^{i+1}} K_{AV}ⁱ (n/nᵢ)ⁱ]
           · (∇C_{(AV)⁺} + C_{(AV)⁺} (nᵢ/n) ∇(n/nᵢ))
```

**反応速度（質量作用則）:**
```
R_{AI}   = [Σᵢ k_{AI}ⁱ K_Iⁱ (n/nᵢ)ⁱ] · (C_{A⁺}C_{I⁰} − C_{(AI)⁺}/K_{A+I}⁰)

R_{AV}   = [Σᵢ k_{AV}ⁱ K_Vⁱ (n/nᵢ)ⁱ] · (C_{A⁺}C_{V⁰} − C_{(AV)⁺}/K_{A+V}⁰)

R_{I+V}  = [Σᵢⱼ k_{I+V}^{ij} K_Iⁱ K_Vʲ (n/nᵢ)^{i+j}]
           · (C_{I⁰}C_{V⁰} − C_{I⁰}* C_{V⁰}*)

R_{AI+V} = [Σᵢⱼ k_{AI+V}^{ij} K_{AI}ⁱ K_Vʲ (n/nᵢ)^{i+j}]
           · (C_{(AI)⁺}C_{V⁰} − K_{A+I}⁰ C_{I⁰}* C_{V⁰}* C_{A⁺})

R_{AI+AV}= [Σᵢⱼ k_{AI+AV}^{ij} K_{AI}ⁱ K_{AV}ʲ (n/nᵢ)^{i+j}]
           · (C_{(AI)⁺}C_{(AV)⁺} − K_{A+I}⁰ K_{A+V}⁰ C_{I⁰}* C_{V⁰}*(C_{A⁺})²)
```

**変数・記号の定義:**

| 記号 | 意味 |
|---|---|
| C_{A⁺}, C_I, C_V, C_{AI}, C_{AV} | 各種の全濃度 |
| C_{I⁰}*, C_{V⁰}* | 点欠陥の熱平衡濃度 |
| K_X | 荷電状態の平衡定数 |
| k_X | 前方反応速度係数 |
| D_X | 各荷電状態の拡散係数 |
| n, nᵢ | 電子濃度, 真性キャリア濃度 |
| i, j | 荷電状態インデックス (−, 0, +) |

---

### 2.3 熱酸化 — Deal-Grove モデル

#### 基礎微分方程式

```
dx₀/dt = B / (A + 2x₀)
```

#### 解析解

```
x₀(t) = (A/2) [√(1 + (t + τ)/τ₀) − 1]

τ₀ = A² / (4B)    (時定数)
τ  = (x_i² + A x_i) / B    (初期酸化膜 x_i による時間オフセット)
```

#### 速度定数の定義

```
B   = 2D C* / N              (放物速度定数 [μm²/hr])
B/A = C* kₛ / N             (線形速度定数  [μm/hr])
A   = 2D (1/kₛ + 1/h)

D  : SiO₂中の酸化剤拡散係数
C* : 酸化剤の SiO₂中平衡溶解度
N  : SiO₂ 1 cm³ あたりの酸化剤分子数（dry O₂: 2.2×10²², wet H₂O: 4.4×10²²）
kₛ : 界面反応速度定数
h  : ガス境界層の物質移動係数
```

#### Arrhenius 温度依存性

```
B     = C₁ exp(−E₁ / kT)
B / A = C₂ exp(−E₂ / kT)
```

| 雰囲気 | C₁ [μm²/hr] | E₁ [eV] | C₂ [μm/hr] | E₂ [eV] |
|---|---|---|---|---|
| 乾式 O₂ (111) | 7.72 × 10² | 1.23 | 6.23 × 10⁶ | 2.00 |
| 湿式 H₂O (111) | 3.86 × 10² | 0.78 | 1.63 × 10⁸ | 2.05 |

#### 圧力依存性

```
Wet:  B(P) = B(1 atm)·P,    B/A(P) = B/A(1 atm)·P
Dry:  B(P) = B(1 atm)·P,    B/A(P) = B/A(1 atm)·P^0.75
```

#### 結晶方位依存性（線形速度定数のみ）

```
(B/A)_{<111>} = 1.68 × (B/A)_{<100>}
```

#### 酸化増速拡散（OED）
熱酸化により Si/SiO₂ 界面で過剰格子間原子が注入される:
```
∂C_I/∂t = ∇(D_I ∇C_I) + G_I − R_I
G_I : 酸化速度に比例した生成項  G_I ∝ dx₀/dt
```

これが後段拡散での B/P 拡散増速（OED）の原因。

---

### 2.4 欠陥・ダメージ・TED

#### Kinchin-Pease 変位数

イオンが反跳原子に転移するエネルギー T に対する変位数 ν(T):
```
ν(T) = 0              (T < Ed)
     = 1              (Ed ≤ T < 2Ed/0.8)
     = 0.8 T / (2Ed)  (T ≥ 2Ed/0.8)
```

Si の変位閾エネルギー: `Ed ≈ 15 eV`

#### NRT（Norgett-Robinson-Torrens）改良版

```
ν_NRT(T) = 0.8 T_dam / (2 Ed)
T_dam = T · (1 − ξ_e)      (電子的エネルギー損失補正)
```

#### "+1" モデル（TED 初期条件）

BCA 各カスケードでフレンケル対が生成されるが、最終的に格子間原子 1 個が残留:
```
ΔC_I(per ion) = Φ / (V_cell × N_ions)     (ドーズあたり 1 格子間原子)
```

この過剰 C_I が TED の駆動力になる。

#### {311} 欠陥・SMIC モデル（Sentaurus 独自）

格子間原子クラスターの動力学（Ostwald 粗大化）:
```
∂C_{311}/∂t = G_{311}(C_I) − E_{311}(C_{311}, T)

C_I^{release} = C_{311} / τ_{311},    τ_{311} = τ₀ exp(Eₐ / kT)
```

この放出された C_I が引き続き TED として拡散。

---

### 2.5 ドーパント活性化・クラスタ

`pdbSet Silicon <Species> ActiveModel <model>` で選択。

| モデル名 | 説明 |
|---|---|
| `None` | 全濃度が即座に活性（デフォルトでない） |
| `Solid` | 固溶度 C_ss(T) でクランプ、超過分は不活性 |
| `Transient` | 時定数付き活性化 |
| `Precipitation` | 析出相の核生成・成長（Ostwald 粗大化） |
| `Cluster` | ドーパント-欠陥クラスタ |
| `BIC` | ボロン格子間クラスタ（B_nI_m） |
| `ChargedCluster` | 荷電状態込みクラスタ |
| `ComplexCluster` | 多成分複合クラスタ |
| `Equilibrium` | 常に化学平衡 |

**BIC（ボロン格子間クラスタ）反応例:**
```
B_n + I ⇌ B_n+1 (または B_nI_m クラスタ形成)
```

**固溶度クランプ（Solid モデル）:**
```
C_active = min(C_total, C_ss(T))
C_ss(T) = C₀ exp(−Eₛₛ / kT)
```

---

### 2.6 応力モデル

粘弾性構成則（Maxwell 体）:
```
dσ/dt = E (dε/dt − dε_th/dt − dε_creep/dt)

dε_creep/dt = σ / (2η)    (クリープ速度、η は粘度)

η(T) = η₀ exp(Eη / kT)    (SiO₂ 粘度の Arrhenius)
```

応力と拡散の結合（応力場による拡散係数変調）:
```
D_eff = D₀ exp(−(Eₐ − Ω·σ_hydro) / kT)    (Ω: 活性化体積)
```

---

### 2.7 境界条件

**HomNeumann（無フラックス、既定）:**
```
J · n̂ = 0
```

**Dirichlet（濃度固定）:**
```
C = C_surface(T)    （例: 固溶度で規定）
```

**偏析（Segregation、界面）:**
```
C₁ / C₂ = m_seg
m_seg = [Arr 0.75 −0.336]  = 0.75 exp(0.336/kT)  （As の SiO₂/Si 界面例）
```

**界面再結合（K_surf モデル）:**
```
J · n̂ = K_surf · (C − C*)
K_surf = A_surf · exp(−E_surf / kT)
```

**Native Oxide デフォルト:** 1.5 nm（全 Si 界面に自動付与）

---

## 3. 数値計算手法

---

### 3.1 空間離散化 — ボックス法（有限体積法）

各制御体積 Ωᵢ で積分形式の連続の式を適用:

```
∂/∂t ∫_{Ωᵢ} C dV + ∮_{∂Ωᵢ} J · n̂ dS = ∫_{Ωᵢ} R dV
```

面 f_{ij} における拡散フラックスの離散化（調和平均拡散係数 + 非直交補正）:
```
∫_{f_ij} J · n̂ dS ≈ D_h · |f_ij| / d_ij · (C_j − C_i) + D_h · (g_f · k_f)

D_h = harmonic_mean(D_i, D_j)    (界面拡散係数)
d_ij = |x_j − x_i|              (セル中心間距離)
k_f  = g_f − (g_f · n̂)n̂        (非直交補正ベクトル)
```

**メッシュ品質指標（直交性）:**
```
orth = (x_j − x_i) · n̂_{ij} / |x_j − x_i|   ∈ (0, 1]
```

---

### 3.2 時間積分

**後退オイラー（1 次、A-安定）:**
```
|Ωᵢ| (C_i^{n+1} − C_i^n) / Δt = F^{n+1}(C^{n+1})
```

**TR-BDF2（2 次、L-安定、大時間刻みに有利）:**
```
段階 1 (γ = 2−√2):
  C^{n+γ} = C^n + γΔt/2 · (F^n + F^{n+γ})
段階 2:
  C^{n+1} = C^{n+γ} + (1−γ)Δt/(2−γ) · (F^{n+γ} + F^{n+1})
```

時間刻み Δt は収束状況に応じて自動制御（エラー推定器ベース）。

---

### 3.3 非線形ソルバー — Newton-Raphson

残差ベクトル **F**(u)=0 に対して:
```
J(u^k) Δu^k = −F(u^k)
u^{k+1}     = u^k + Δu^k

J = ∂F/∂u   (ヤコビアン行列、解析的に構築)
```

収束判定:
```
||Δu^k||_∞ / (||u^k||_∞ + ε) < tol_Newton
```

**Five-stream では 5 変数（C_A, C_I, C_V, C_AI, C_AV）を完全結合した
ブロック構造のヤコビアンを用いる（Picard では収束しない）。**

---

### 3.4 線形ソルバー

**直接法（大規模 2D/3D メイン）:**
- PARDISO（Intel MKL）: LU 分解、ネスティド-ダイセクション ordering
- UMFPACK: 研究用代替

**反復法（補助・大規模 3D）:**
- ILU(0)/ILUT 前処理 GMRES
- AMG（代数的多重格子）

---

### 3.5 Arrhenius パラメータ評価

Sentaurus Process は実行時に温度 T で評価:

```tcl
[Arrhenius A Ea]          → A · exp(−Ea / (kT))   [eV 単位]
[Arr A Ea]                → 同上（短縮形）
[pdbDelayDouble Reg Sp Par]  → 未評価式を取得
[pdbGet Reg Sp Par]          → 現在温度で評価した値を返す
```

パラメータの上書き:
```tcl
pdbSet Silicon Boron Dstar {[expr 2.0 * [pdbDelayDouble Silicon Boron Dstar]]}
pdbSet Oxide_Silicon Vacancy Ksurf {[Arrhenius 1e3 0.1]}
```

---

## 4. パラメータDB

### 階層構造

```
Region / Interface
  └─ Species
       └─ Parameter = Value|Arrhenius 式
```

例:
```tcl
pdbSet  Silicon        Boron      Dstar         {[Arr 0.037 3.46]}
pdbSet  Silicon        Boron      ActiveModel   Solid
pdbSet  Oxide_Silicon  Arsenic    Segr          {[Arr 0.75 -0.336]}
pdbSet  Oxide_Silicon  Vacancy    Ksurf         {[Arrhenius 1e3 0.1]}
pdbSet  Oxide          Grid       perp.add.dist  1e-7
pdbSetBoolean  Grid    DoNotMove.Reaction        1
pdbSet  InfoDefault    <0|1|2>    # 出力冗長度
```

パラメータファイルの格納場所（インストール環境）:
```
$STROOT/tcad/$STRELEASE/lib/score/Params/
```

> **Cprocess 方針:** pdbSet による実行時上書きは再現しない。代わりに
> `materials.cpp` の係数を変更するか、将来的に外部 YAML/JSON で管理する。

---

## 5. 単位系・座標系

### 単位（角括弧 `<>` でコマンド内に明示）

| 物理量 | Sprocess 単位 | 角括弧省略時の既定 |
|---|---|---|
| 長さ | `<um>`, `<nm>`, `<cm>` | µm |
| エネルギー | `<keV>`, `<eV>`, `<MeV>` | keV |
| 濃度 | `<cm-3>` | — |
| ドーズ | `<cm-2>` | — |
| 時間 | `<min>`, `<s>`, `<hr>` | min (diffuse) |
| 温度 | `<C>`, `<K>` | °C |
| 圧力 | atm | 1 atm |

### 座標系

```
Sprocess (1D/2D):  x  ↓ (ウェハ深さ、下向き正)
                   y  → (横方向)
Cprocess:          z  ↓ (ウェハ深さ)
                   x, y (横方向)
```

変換ルール: Sprocess の x が Cprocess の z に対応。

---

## 6. 参照

- Sentaurus Process チュートリアル 1D (sp_2)
  <https://ghzphy.github.io/Sentaurus_Training/sp/sp_2.html>
- Sentaurus Process チュートリアル 2D (sp_3)
  <https://ghzphy.github.io/Sentaurus_Training/sp/sp_3.html>
- Sentaurus Process チュートリアル モデル/pdbSet (sp_4)
  <https://ghzphy.github.io/Sentaurus_Training/sp/sp_4.html>
- Sentaurus Process チュートリアル SiC (sp_15)
  <https://ghzphy.github.io/Sentaurus_Training/sp/sp_15.html>
- Five-Stream Dunham Model（TU Wien Hollauer 博士論文 §4.2）
  <https://www.iue.tuwien.ac.at/phd/hollauer/node24.html>
- Dual-Pearson 分布（TU Wien Khalil 博士論文 §5.3.1）
  <https://www.iue.tuwien.ac.at/phd/khalil/node64.html>
- Deal-Grove 酸化モデル（TU Wien Hollauer §2.6）
  <https://www.iue.tuwien.ac.at/phd/hollauer/node16.html>
- B.E. Deal & A.S. Grove, J. Appl. Phys. 36 (1965) 3770
- R.B. Fair, "Impurity Doping Processes in Silicon" (1981) — 拡散係数
- J.F. Ziegler, J.P. Biersack, U. Littmark, "The Stopping and Range of
  Ions in Solids" (1985) — ZBL ポテンシャル
- Norgett-Robinson-Torrens, Nucl. Eng. Design 33 (1975) — NRT モデル
- Dunham & Wu, J. Appl. Phys. 78 (1995) — Five-Stream モデル
