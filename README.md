# Cprocess

3次元半導体プロセスシミュレータ(C++17、外部依存なし)。

非構造四面体メッシュ上のセル中心**有限体積法 (FVM)** で、イオン注入と不純物拡散を解きます。Gmsh メッシュの読み込みと ParaView 用 VTU 出力に対応しています。

## 特徴

- **非構造四面体メッシュ**
  - Gmsh ASCII `.msh`(フォーマット 2.2 / 4.1、1次要素)読み込み。Physical Volume → 領域(材料)、Physical Surface → 境界パッチ
  - 内蔵ボックスメッシャ(六面体を整合的な Kuhn 6分割で四面体化)
- **FVM 離散化**
  - 2点フラックス近似(TPFA)+ 面拡散係数の調和平均
  - 重み付き最小二乗 (LSQ) セル勾配による**非直交補正**(遅延補正方式)。LSQ 勾配は任意の四面体メッシュ上で線形場に対し厳密
  - 後退オイラー時間積分(無条件安定)+ 濃度依存拡散係数に対する Picard 反復
  - Jacobi 前処理付き CG(対称系)、フォールバックに BiCGSTAB
- **プロセス工程**
  - イオン注入(解析): ガウス分布(深さ方向)+ マスク窓の誤差関数畳み込み(横方向広がり)。B/P/As/Sb の近似レンジテーブル内蔵、`rp=`/`drp=` で直接指定も可
  - イオン注入(モンテカルロ / BCA): TRIM 方式の二体衝突近似。ZBL 万有ポテンシャルの**厳密散乱積分**を事前テーブル化して参照(corteo 方式)、Lindhard–Scharff 電子的阻止能、チルト/ローテーション対応。**スレッド並列**かつチャンク分割 RNG +決定的リダクションにより**スレッド数によらずビット同一の結果**
  - 拡散: Fair の荷電点欠陥モデル(温度・キャリア濃度依存 D)、電界増速効果、複数不純物の n/ni を介した連成
  - 析出(predeposition): 境界パッチへの固定表面濃度(Dirichlet)境界条件
  - 活性化濃度: 保存時に固溶度でクランプした `*_active` と `NetDoping` を出力
- **入出力**: SUPREM 風のテキストデッキ入力、VTK XML (`.vtu`) 出力

## ビルド

```bash
cmake -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j
ctest --test-dir build        # 検証テスト実行
```

要件: C++17 コンパイラと CMake 3.16+ のみ。

## クイックスタート

```bash
./build/cprocess examples/pwell.in
```

`pwell.vtu` が生成されます。ParaView で開き、`NetDoping` をセルデータとして表示すると p-well 接合が確認できます(座標は µm 単位で出力されます)。

```
# examples/pwell.in の中身(抜粋)
mesh box xmax=1um ymax=1um zmax=1.5um nx=12 ny=12 nz=40
init species=P conc=1e15                  # n型基板
implant species=B energy=60keV dose=2e13 x1=0.25um x2=0.75um y1=0.25um y2=0.75um
diffuse time=60min temp=1100C dt=60s      # ドライブイン
save file=pwell.vtu
```

## デッキコマンド

1行1コマンド、`#` 以降はコメント。引数は `key=value` 形式。

| コマンド | 説明 |
|---|---|
| `mesh box xmax= ymax= zmax= nx= ny= nz= [xmin= ymin= zmin=]` | ボックスメッシュ生成。パッチ名は `xmin,xmax,ymin,ymax,zmin,zmax` |
| `mesh gmsh file=dev.msh [scale=1um]` | Gmsh メッシュ読み込み。`scale` はメッシュ座標 1 単位の実寸 |
| `region {tag=N \| name=NAME \| all} material=silicon\|oxide\|nitride\|poly\|gas` | 領域への材料割り当て(シリコン以外は注入・拡散の対象外の不活性領域) |
| `init species=B conc=1e15 [region=N]` | 一様背景濃度の設定 |
| `implant species=P dose=1e13 {energy=80keV \| rp=0.1um drp=0.04um} [drl=] [x1= x2= y1= y2=]` | ガウス注入。`x1..y2` を全部与えるとマスク窓、`drl` は横方向標準偏差(既定 0.8ΔRp) |
| `implant species=B dose=1e13 energy=50keV method=mc [ions=100000] [threads=0] [seed=1] [tilt=7] [rotation=30] [x1=..y2=]` | モンテカルロ注入。`threads=0` で全コア使用。同一 seed なら結果はスレッド数に依存しない |
| `bc species=P patch=zmax conc=1e20` / `bc clear` | 表面濃度固定(析出)境界条件の設定/解除。以後の `diffuse` に適用 |
| `diffuse time=30min temp=1000C [dt=30s] [fieldenh=on\|off] [nonortho=on\|off]` | 拡散アニール。存在する全不純物を連成して解く |
| `save [file=out.vtu]` | VTU 出力(各不純物、`*_active`、`NetDoping`、`Region`) |
| `print` | 各不純物の積分量・ピーク濃度などを表示 |
| `stop` | デッキ実行終了 |

**単位**: 接尾辞 `um/nm/mm/cm/m`(長さ、既定 cm)、`s/min/h`(時間)、`C/K`(温度、既定 K)、`keV/eV/MeV`(エネルギー)。濃度は cm⁻³、ドーズは cm⁻² の無接尾辞数値。

## 物理モデル

### イオン注入(解析ガウス)

垂直入射、平坦上面(メッシュ z 最大面)を仮定:

$$C(d) = \frac{\Phi}{\sqrt{2\pi}\,\Delta R_p} \exp\!\left(-\frac{(d-R_p)^2}{2\Delta R_p^2}\right), \quad d = z_{top} - z$$

マスク窓 $[x_1,x_2]\times[y_1,y_2]$ がある場合は横方向に誤差関数で畳み込み:

$$C \mathrel{*}= \tfrac{1}{2}\!\left[\mathrm{erf}\tfrac{x-x_1}{\sqrt{2}\Delta R_l} - \mathrm{erf}\tfrac{x-x_2}{\sqrt{2}\Delta R_l}\right] \cdot (\text{y も同様})$$

内蔵の $R_p/\Delta R_p$ テーブルは LSS 理論ベースの教科書近似値です。精度が必要な場合は SRIM 等の値を `rp=`/`drp=` で直接指定してください。

### イオン注入(モンテカルロ / BCA)

`method=mc` で TRIM 方式の二体衝突近似(BCA)によるイオン履歴追跡を行います(アモルファスターゲット):

- **自由飛行**: 原子間距離 $L = N^{-1/3}$ ごとに 1 回衝突。衝突径数は $p = p_{max}\sqrt{U}$($\pi p_{max}^2 L N = 1$、$U$ は一様乱数)
- **核的散乱**: ZBL 万有遮蔽ポテンシャルの**厳密古典散乱積分**

  $$\theta = \pi - 2\int_{R_0}^{\infty} \frac{b\,dR}{R^2\sqrt{1 - V(R)/\varepsilon - b^2/R^2}}$$

  を起動時に $(\log\varepsilon, \log b)$ の 160×160 テーブルへ事前計算し、双線形補間で参照します(corteo 方式)。換算単位では全イオン種・ターゲットで同一関数なのでテーブルは 1 個。MAGIC 公式のようなフィット式を使わないため、フィット係数由来の誤差がありません(テストで独立実装の数値積分と < 2% 一致を検証)
- **電子的阻止能**: Lindhard–Scharff $S_e = k\sqrt{E}$($k$ は $Z_1, Z_2, M_1$ から決定)。≲ 数百 keV で妥当、500 keV 超で警告
- **エネルギー移行**: $T = T_{max}\sin^2(\theta_c/2)$、$T_{max} = 4M_1M_2E/(M_1+M_2)^2$。実験室系偏向角は質量比から厳密変換
- **停止**: $E < 5$ eV で静止とみなし、静止点を含む四面体セル(背景格子による点位置検索)へドーズ重み $\Phi A_{source}/N_{ions}$ で計上
- **境界**: 上面流出=後方散乱、下面流出=透過として計数。全面注入では横方向に周期境界、窓注入では領域外流出を計数
- **チルト/ローテーション**: ビーム方向を `tilt=`/`rotation=`(度)で指定可

**並列設計**: イオン履歴を 2048 個単位のチャンクに分割し、チャンクごとに counter-based シードの独立 RNG ストリーム(splitmix64 → mt19937_64、一様乱数はビット演算で生成しライブラリ実装非依存)を割り当てます。ワーカーはスレッドローカルな**整数**ヒットカウントとチャンク別モーメント部分和に蓄積し、最後に逐次リダクション — 浮動小数の加算順序が固定されるため、**結果はスレッド数・スケジューリングによらずビット同一**です。履歴ループは共有状態を持たず(散乱テーブル・セル探索器は構築後 immutable)、将来の MPI / GPU 展開にもそのまま載る構造です。

### 不純物拡散(Fair モデル)

$$\frac{\partial C}{\partial t} = \nabla\!\cdot\!\left(h\, D(T, n/n_i)\, \nabla C\right)$$

$$D = D^0 + D^-\frac{n}{n_i} + D^=\left(\frac{n}{n_i}\right)^2 + D^+\frac{p}{n_i}, \qquad D^X = d_X e^{-E_X/kT}$$

- 係数は Fair の総説(1981)の標準値(B, P, As, Sb)
- $n/n_i$ は全不純物の正味ドーピングから電荷中性で決定(複数不純物が D を介して連成)
- $n_i(T) = 3.87\times10^{16}\, T^{1.5} e^{-0.605\,\mathrm{eV}/kT}$ cm⁻³(Morin–Maita)
- 電界増速係数 $h = 1 + \frac{|c|}{\sqrt{1+c^2}} \le 2$($c = N_{net}/2n_i$、支配的キャリアと同型の不純物に適用)。`fieldenh=off` で無効化
- 固溶度(近似 Arrhenius フィット)は出力時の活性濃度クランプにのみ使用

### FVM 離散化

面 $f$(セル P–N 間)のフラックス:

$$F_f = \underbrace{D_f\, \frac{|\mathbf{S}|^2}{\mathbf{S}\cdot\mathbf{d}}\,(C_N - C_P)}_{\text{陰的 (TPFA)}} + \underbrace{D_f\, \nabla C_f \cdot \mathbf{k}}_{\text{陽的非直交補正}}, \quad \mathbf{k} = \mathbf{S} - \frac{|\mathbf{S}|^2}{\mathbf{S}\cdot\mathbf{d}}\mathbf{d}$$

- $D_f$: 面法線方向距離による調和平均(材料界面・急峻な D 勾配で頑健)
- $\nabla C_f$: 重み付き最小二乗セル勾配の面内挿。LSQ は線形場で厳密なため、補正後のスキームは任意の四面体形状で整合(下記検証参照)
- 境界閉包: Dirichlet 面はゴースト点、ゼロフラックス面は $\nabla C\cdot\hat{n}=0$ 制約行として LSQ に組み込み
- 時間積分は後退オイラー。行列は対称 M 行列となり CG で解く

## 検証(`ctest`)

| テスト | 内容 |
|---|---|
| `mesh` | セル体積の総和、面ベクトル閉包 $\sum_f \mathbf{S}_f = 0$、パッチ面積、反転四面体の自動修正 |
| `solver` | CG / BiCGSTAB の製造解収束 |
| `gmsh` | MSH 2.2 / 4.1 の読み込み、物理グループ → 領域/パッチ対応、スケール変換 |
| `diffusion` | (1) 埋め込みガウスの分散成長 $\sigma^2 += 2Dt$(解析解比 5% 以内 — 実測 ~0.06%、アスペクト比 5 の歪んだメッシュ上)、質量保存 10⁻⁶ 以内。(2) 析出の erfc プロファイル(L2 誤差 < 10%)。(3) 材料界面のゼロフラックス(漏れなし) |
| `mc` | (1) 散乱テーブル vs 独立実装の厳密散乱積分(< 2%)。(2) B 50 keV の $R_p/\Delta R_p$ が物理レンジ内・収支の厳密一致(deposited+lost = ions)・ビニングのドーズ/モーメント整合。(3) **スレッド数 1 と 4 でビット同一**。(4) 窓注入の横方向封じ込め |

非直交補正の効果: アスペクト比 10 の Kuhn 四面体メッシュで、補正なし TPFA は分散成長を 33% 過小評価しますが、LSQ 補正ありでは 0.06% に収まります(メッシュ細分化に対して不変)。

## 制限事項と今後の拡張

現バージョンの主な仮定・制限:

- 解析注入は垂直入射のみ。MC 注入はチルト対応だが、上面(z 最大面)が平坦であることを仮定
- MC 注入はアモルファスターゲット(チャネリングなし)。全領域にシリコンの阻止能を適用(酸化膜マスク等の材料別阻止能は未対応)。反跳カスケード・損傷蓄積は追跡しない
- 酸化(Deal–Grove、移動境界)、エッチング/堆積(メッシュ変形)は未実装
- 点欠陥連成(TED)、クラスタリング/活性化動力学、材料界面の偏析は未実装(固溶度クランプは後処理のみ)
- シリコン以外の材料は拡散に関して不活性
- 拡散の線形ソルバは逐次実行(OpenMP / ILU0 / AMG は今後)。MC 注入はスレッド並列済み

## ライセンス

リポジトリ同梱の LICENSE を参照してください。
