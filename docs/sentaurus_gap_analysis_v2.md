# Cprocess 総合評価と Sentaurus Process 残存ギャップ分析 (v2)

作成日: 2026-07-14 / 対象: HEAD `308c2ae`(全ロードマップ P1-1〜11, P2-1〜8, S-1〜5, M-1〜6, PA-1〜5, P3-a〜h 実装完了時点)

本書は初版ギャップ分析(`sentaurus_gap_analysis.md`、実装前に作成)の後継である。初版で「大」とされた
ギャップ(点欠陥/TED、酸化+移動境界、活性化、応力、ParamDB、構造編集)は**全て何らかの形で実装済み**
となったため、本書では (1) 現時点の到達度の総合評価、(2) なお残る SProcess との差分の構造的な整理、
(3) 次に投資すべき領域の優先順位付け、を行う。評価は 3 系統の独立監査
(物理モデル / 数値・メッシュ・並列基盤 / フロントエンド・I/O・校正)に基づき、
すべての判定は実コード(file:function)を一次情報として検証した。

凡例: ✓ 実装済み / △ 簡略化・部分実装・未統合 / ✗ 未実装。
深刻度は実用 TCAD(ロジック/パワー系プロセス設計)観点での影響度。

---

## 0. 総合評価(エグゼクティブサマリ)

**到達点**: 注入(解析 + 化合物ターゲット MC-BCA + チャネリング + ダメージ)→ SPER → RTA/TED
(3 場点欠陥 + {311} + クラスタ)→ 酸化(1D Deal-Grove + 2D LOCOS + OED + 応力遅延)→
エピ/シリサイド/エッチ・デポ(レベルセット)→ 力学(FEM + SiGe 歪み)→ デバイス連携出力、
という**単一パイプラインが物理的に連成して一気通貫で回る**。テストは C++ 49 本 + Python 51 本が
全て通り、3 層アーキテクチャ(C++ core / proc:: / Python)と単位規約は一貫している。
初版分析時点の「注入と拡散しかない骨格」からは質的に別物である。

**構造的な残差**: 各領域とも SProcess の「校正済み・完全結合」実装に対して 1 段簡略な
「実効モデル」に留まる。残るギャップは次の 4 カテゴリに集約される。

| カテゴリ | 本質 | 代表例 |
|---|---|---|
| A. 校正データの薄さ | モデルの形はあるが係数が教科書級/テスト合わせ | 注入モーメント表(8種×5–8点)、クラスタ Eb、ox2d.nitride_leak |
| B. 結合度の一段の浅さ | 演算子分割/実効場で代替 | 五流(React)不在、荷電欠陥なし、粘性流動酸化なし |
| C. 部品はあるが未統合 | 実装済みコンポーネントが主経路に配線されていない | AMG、BCSR ブロック解法、GPU CG、coarsen、既定 1 スレッド化 |
| D. フロントエンドの遅れ | エンジンに UI/文書/変換器が追いついていない | sprocess.py 変換器の陳腐化、デッキ 16 vs Python 37、旧文書の誤記 |

このうち **C と D は工数小・効果大**(配線と更新だけでよい)、**A は地道なデータ整備**、
**B はアーキテクチャ投資**を要する。

---

## 1. 物理モデル: 到達度と残差

### 1.1 イオン注入

| 機能 | SProcess | Cprocess | 根拠 | 深刻度 |
|---|---|---|---|---|
| Pearson-IV 解析分布 | ✓(校正済み) | ✓ ODE 数値構築、不正モーメント時 Gauss フォールバック | implant.cpp:build_pearson4 | 低 |
| Dual-Pearson(チャネリングテール) | ✓ 9パラメータ + 校正ライブラリ | ✗ 単峰のみ | — | **高** |
| モーメントテーブル | 全種×広エネルギー実測 DB | △ 8種×5–8点(10–200 keV、教科書級、範囲外クランプ) | materials.cpp:kDopants | **高** |
| MC-BCA(ZBL/化合物/チャネリング/KP ダメージ/動的非晶質化) | ✓ Crystal-TRIM | ✓ P3-a で化合物ターゲット(SiO2/Si3N4 成分別衝突)まで到達 | mc_implant.cpp | 低 |
| 注入間の累積非晶質化 | ✓ | △ damage 場は SPER 用に永続化されるが、次の MC 注入の結晶状態には**戻らない**(毎回フレッシュ) | process.cpp:implant_mc | 中 |
| 分子注入(BF2 等) | ✓ | ✗ | — | 中 |
| 一般段差の遮蔽(shadowing) | ✓ | △ レジストスタック透過のみ | mc_implant.cpp | 中 |

### 1.2 拡散・アニール・活性化

| 機能 | SProcess | Cprocess | 根拠 | 深刻度 |
|---|---|---|---|---|
| Fermi 4 荷電項 D + 電界増速 + 電荷中性 Picard | ✓ | ✓ | materials.cpp:dopant_diffusivity | 低 |
| 点欠陥陽解(I/V) | Pair/ChargedPair/React 階層 | △ C_I・C_V・C_311 の 3 場(中性のみ)。ドーパントは (1−fi)+fi·(C_I/C_I*) の準平衡 Pair 相当 | diffusion.cpp:step_once_ted | 中 |
| 五流(React)/AI・AV ペア方程式/荷電欠陥状態 | ✓ Newton 完全結合 | ✗ | — | **高** |
| TED({311} 捕獲/放出) | ✓ + 転位ループ/EOR | △ 単一実効 {311} 場。EOR 欠陥帯の空間構造なし | diffusion.cpp | 中 |
| クラスタ活性化 | 9 モデル(BIC/Charged/析出/Transient…) | △ B(BIC)/As(As4V)/C(C-I シンク)のみ。P/Sb/In なし。Eb はテスト合わせ校正(コメント明記) | materials.cpp:cluster_params | **高** |
| 固溶度クランプ + SPER 準安定活性 | ✓ | ✓(act_factor×C_ss) | process.cpp:sper | 低 |
| 偏析・ドーズロス(多材料 SegTable + シリサイド) | ✓ | ✓ | diffusion.hpp:SegTable | 低 |
| RTA ランプ / 適応 dt | ✓ | ✓(区分線形 / step-doubling opt-in) | DiffuseOpts | 低 |

### 1.3 酸化

| 機能 | SProcess | Cprocess | 根拠 | 深刻度 |
|---|---|---|---|---|
| Deal-Grove 1D(dry/wet、44/56% 体積収支、メッシュ移動) | ✓ | ✓ | oxidation.cpp, process.cpp:oxidize | 低 |
| Massoud 薄膜補正 | ✓ | ✗(native 25Å オフセットのみ) | — | 中 |
| 圧力/HCl/方位依存 | ✓ | ✗(1 atm、<100> 固定) | — | 中 |
| 2D LOCOS/バーズビーク | 粘弾性流動 + 移動境界 ALE | △ 酸化剤の定常拡散-反応(Robin)は本物だが、幾何は (x,y) カラム毎 1D 界面高さ + `ox2d.nitride_leak` チューニング係数。粘性流動なし、box メッシュ限定 | process.cpp:oxidize_2d | **高** |
| OED | ✓ | ✓(界面 Si 消費フラックス比例 I 注入) | process.cpp:oxidize | 低 |
| 応力遅延酸化 | ✓ | △ P3-e opt-in(ks の指数スケール) | oxidant_solver.hpp | 中 |
| poly/SiGe 酸化 | ✓ | ✗(bare Si のみ、他は throw) | — | 中 |

### 1.4 力学・応力・SiGe

| 機能 | SProcess | Cprocess | 根拠 | 深刻度 |
|---|---|---|---|---|
| FEM(P1 tet 線形弾性 + 固有ひずみ: 熱/膜真性/SiGe Vegard) | 粘弾性 FEM | ✓(弾性)/△(粘弾性は 1 ステップ Maxwell 減衰のみ) | fem.cpp, mechanics.cpp | 中 |
| 応力→拡散(exp(−pV_act/kT)) | ✓ | ✓ P3-e(B/P/As の V_act、opt-in) | diffusion.hpp:pressure | 低 |
| 応力履歴・境界移動との連成 | ✓ | △ mechanics はユーザ明示呼出しの一発解 | process.cpp:mechanics | 中 |
| Si 異方性弾性 | ✓ | ✗(等方のみ) | — | 低 |

### 1.5 シリサイド・エピ・材料

| 機能 | SProcess | Cprocess | 根拠 | 深刻度 |
|---|---|---|---|---|
| シリサイド | 多相系列(Ni2Si→NiSi)、2D | △ 単相 NiSi/TiSi2、x²=x0²+Bt、ブランケットのみ、偏析ドーズロスあり | process.cpp:silicide | 中 |
| エピタキシー | 成長速度モデル、選択エピ、facet | △ 幾何成長 + in-situ 一様ドープ + 自動熱履歴 | process.cpp:epitaxy | 中 |
| ドーパント種 | B,P,As,Sb,In,Ga,Al,C,F,Ge,N,H + 分子 | △ 8 種(B,P,As,Sb,In,C,F,Ge) | materials.cpp:kDopants | 低〜中 |
| 基板材料 | Si,SiO2,Si3N4,poly,SiGe,SiC(4H/6H/3C),Ge,III-V | △ Si,SiO2,Si3N4,poly,NiSi/TiSi2,Ni/Ti,gas,resist。**SiC/Ge 基板なし**(SiGe は Si 中 Ge 場の歪み源のみ) | materials.hpp:MatId | 中(SiC 用途では致命的) |

---

## 2. 数値・メッシュ・並列基盤: 到達度と残差

### 2.1 離散化・ソルバー

| 機能 | SProcess | Cprocess | 根拠 | 深刻度 |
|---|---|---|---|---|
| 制御体積 | ボックス法(Voronoi 双対) | △ セル中心 FVM + 二点フラックス + 遅延非直交補正 + M 行列クランプ | diffusion.cpp:build | 中 |
| Kuhn-tet 対角面の非直交問題 | 該当なし | △ **既知の構造的問題**。定常 Laplace では遅延補正が不安定 → oxidant ソルバーは hex 集約(box メッシュ専用)で回避。一般非構造メッシュの定常拡散は解けない | oxidant_solver.cpp 冒頭 | **高** |
| 疎直接法(PARDISO 級) | ✓ 主力 | ✗(ソルバー失敗が即 throw になる脆弱性の根源) | — | **高** |
| Krylov(CG/BiCGSTAB/GMRES)+ ILU(0) レベル並列 | ✓ | ✓ 主経路に統合済み | sparse.cpp | 低 |
| AMG | ✓ | △ 2 レベル SA 実装済みだが**主経路未統合**(テストのみ) | sparse.cpp:cg_amg | 中 |
| ブロック解法(BCSR + block Jacobi) | 完全結合の基盤 | △ 実装済みだが **run_ted が未使用**(逐次分離のまま) | sparse.cpp:bicgstab_bjacobi | **高** |
| RCM リオーダリング | (ND in PARDISO) | ✓ SpMV 局所性目的で統合済み | sparse.cpp:rcm_order | 低 |
| 完全結合多種 Newton | ✓ 解析ヤコビアン(five-stream 必須要件) | ✗ JFNK(S-3)は種内アフィン部分問題の実質線形解。種間は Picard | diffusion.cpp コメント | **高** |
| 時間積分 | BE + TR-BDF2 | △ BE のみ + step-doubling 適応 dt(opt-in) | — | 中 |

### 2.2 メッシュ

| 機能 | Sentaurus Mesh | Cprocess | 根拠 | 深刻度 |
|---|---|---|---|---|
| 初期メッシュ | 境界適合 Delaunay | △ box + Kuhn 6-tet、Gmsh 取込 | box_mesh.cpp | 中 |
| プロセス中の自動適応リメッシュ | ✓ refinebox | △ refine(勾配/異方性)は**手動・操作単位**。coarsen は実装済みだが **proc::/Python 未配線** | remesh.cpp | **高** |
| 品質修復(平滑化/分割/フリップ) | メッシュ側で保証 | ✓ ALE 移動後に自動 | remesh.cpp:repair_quality | 低 |
| 界面適合 | 常時境界適合 | △ 消極的保護のみ(界面=region タグの階段状) | remesh.cpp | 中 |
| レベルセット | ✓ | △ 1 次アップウィンド + 再初期化(低次) | levelset.cpp | 中 |
| 汎用フィールド転写 | 保存的補間 | △ 最近傍 O(N²)・非保存 | field_transfer.cpp | 中 |
| ネイティブ 1D/2D | ✓ | △ 退化 3D(mesh1d) | box_mesh.cpp | 低〜中 |

### 2.3 並列

| 機能 | SProcess | Cprocess | 根拠 | 深刻度 |
|---|---|---|---|---|
| OpenMP(SpMV/ILU レベル/面彩色アセンブリ/種並列) | ✓ | ✓(決定性設計込み) | sparse.cpp, diffusion.cpp | 低 |
| **既定スレッド数の抑制** | — | ⚠ OMP_NUM_THREADS 未設定時に既定 1 スレッド + run_ted 全体シングルスレッド化(CI 環境対策)。**実マルチコアの既定性能を殺している** | sparse.cpp:ensure_sane_thread_count | 中 |
| GPU | — | △ cg_jacobi_gpu 足場のみ(拡散未接続) | CPROCESS_GPU | 中 |
| MPI | ✓ | △ MC チャンク分散のみ | CPROCESS_MPI | 中 |

---

## 3. フロントエンド・I/O・校正: 到達度と残差

### 3.1 コマンド/スクリプト面

- テキストデッキ: **16 コマンド**(mesh/region/init/implant/photo/mask/strip/bc/diffuse/oxidize/epitaxy/deposit/silicide/save/print/stop)
- Python `Simulation`: **37 メソッド**(+pybind 44)— etch/etch_rate/deposit_conformal/oxidize_2d/sper/mechanics/refine/mask_polygon/load_gds/ramp/set_param/save_state/export_device 等は **Python 専用でデッキから呼べない**
- SProcess Tcl の制御構文(変数/if/foreach)相当はデッキになし(Python が代替)

| 機能 | SProcess | Cprocess | 深刻度 |
|---|---|---|---|
| init tdr=(TDR 読込)/TDR 出力 | ✓ | ✗(独自 CPRC1 + VTU + meta.json) | **高** |
| contact 定義 → SDevice 直結 | ✓ | ✗(export_device は境界パッチ名のみ) | **高** |
| 抽出(Xj/Rs/layers/WritePlx/select DSL) | ✓ | △ dose()/profile()/print のみ | 中 |
| 名前付きマスク/negative/segments | ✓ | △ 矩形窓 + mask_polygon + **GDSII 読込(独自の強み)** | 中 |
| temp_ramp/雰囲気(分圧/gas_flow)のデッキ対応 | ✓ | △ Python ramp= のみ、雰囲気なし | 中 |
| transform reflect/rotate/cut | ✓ | ✗ | 中 |
| wafer.orient/miscut/SiC | ✓ | ✗ | 低〜中 |

### 3.2 校正(ParamDB)

- **チューナブル**: 拡散 4 荷電項、固溶度、fi、偏析、多材料 D、シリサイド成長則、点欠陥、TED、SPER、OED、2D 酸化、力学、応力結合、SiGe — 主要物理は概ね実行時上書き可能(初版ギャップ「大」は解消)
- **ハードコードのまま**: Deal-Grove B/A、注入モーメント表、MC 阻止能/ZBL テーブル、ni(T)
- **フロー未整備**: 階層 PDB なし(フラット KV)、校正ファイル一括ロードなし、デッキから設定不可、**未知キーが無警告で不活性**(誤字検出不能 — 校正作業で危険)

### 3.3 変換器・文書の陳腐化(即時対応可能な最大の落差)

- `python/cprocess/sprocess.py`(SProcess 入力変換器)が 2026-06-15 時点のエンジン能力を前提:
  **oxidize/deposit/etch/silicide/photo/strip/temp_ramp/pdbSet を「未対応」と切り捨てる**が、実際は全てエンジンにある。対応種も B/P/As/Sb のみ(エンジンは 8 種)。変換器のテストも無い。
- `docs/sprocess_command_map.md` と `docs/sentaurus_gap_analysis.md`(初版)も同様に、実装済み機能を未実装と記載。
- その他: `Simulation.bbox()` が µm API の中で唯一 cm を返す罠、体系的ユーザガイドなし、CI 定義なし。

---

## 4. 残存ギャップ 総合ランキング(全領域統合)

工数と効果で 3 群に分けた総合優先順位。

### 第 1 群: 配線・更新だけで効く(工数小・効果大)

| # | ギャップ | 内容 |
|---|---|---|
| 1 | **sprocess.py 変換器の全面更新** | エンジン実装済みの oxidize/deposit/etch/silicide/photo/strip/ramp/set_param への対訳を追加、8 種対応、変換器テスト新設。SProcess レシピ互換度が最小工数で最大改善 |
| 2 | **陳腐化文書の更新** | sprocess_command_map.md / sentaurus_gap_analysis.md(初版)を本書で置換・注記。ユーザ/後続開発の誤判断を防ぐ |
| 3 | **実装済み部品の主経路への統合** | AMG を DiffusionSolver の前処理選択肢に、BCSR ブロック解法を run_ted の I-V-ドーパント結合に、coarsen を proc::/Python に配線、GPU CG の接続判断 |
| 4 | **既定スレッド数の見直し** | ensure_sane_thread_count の 1 スレッド既定と run_ted の OmpThreadGuard(1) を環境検知 or 明示設定に変更。実機性能の解放 |
| 5 | **デッキと Python の機能同期** | etch/oxidize_2d/sper/mechanics/refine/ramp/set_param 等をデッキコマンド化。ParamDB の未知キー警告追加。bbox() 単位修正 |

### 第 2 群: データ整備(工数中・定量精度に直結)

| # | ギャップ | 内容 |
|---|---|---|
| 6 | **注入モーメント表の拡充 + dual-Pearson** | 解析注入は日常使用の大半を占める。SRIM 等による広エネルギー域(1 keV〜数 MeV)の表整備と双峰分布の導入が実効精度を最も左右する |
| 7 | **クラスタ/活性化パラメータの文献校正** | テスト合わせの Eb(B 2.7/As 2.6 eV)を文献値+検証データで再校正、P/Sb/In へ拡張。SIMS/Rs 実測との突合フローを作る |
| 8 | **薄膜酸化(Massoud)/圧力/HCl/方位依存 + Deal-Grove 係数の ParamDB 化** | ゲート酸化・高圧酸化レシピへの対応。DG 係数のハードコード解消も同時に |

### 第 3 群: アーキテクチャ投資(工数大・SProcess の三本柱)

| # | ギャップ | 内容 |
|---|---|---|
| 9 | **完全結合多種 Newton + 疎直接法** | 五流(React)/荷電欠陥への道。BCSR は在庫あり。解析ヤコビアン + 直接法フォールバック(PARDISO 級 or 外部ライブラリ)が前提 |
| 10 | **2D/3D 酸化の粘性流動 + 一般移動境界** | カラム 1D 界面 + nitride_leak 係数から、粘性流動(Stokes 近似)+ ALE への移行。LOCOS/STI 形状の定量化に必須 |
| 11 | **プロセス中の自動適応リメッシュ + 境界適合** | refine/coarsen/repair の自動オーケストレーション(refinebox 相当)と、階段状界面からの脱却(界面適合 or カットセル) |
| 12 | **Kuhn-tet 非直交性の根本解決** | Voronoi/ボックス法離散化への移行 or 一般メッシュでの安定な非直交補正。hex 集約回避策の box メッシュ依存を解消 |
| 13 | **TDR 互換 or 実用的デバイス連携** | contact 定義 + (可能なら)TDR 読み書き。不可能なら meta.json 経路のリファレンス実装(SDevice 側読み込みスクリプト)整備 |
| 14 | **五流物理そのもの**(荷電 I/V、AI/AV ペア、EOR/転位ループ)、分子注入(BF2)、注入間累積非晶質化の BCA 反映、SiC/Ge 基板 |

---

## 5. 結論

- 初版ギャップ分析の「大」項目は全て解消され、Cprocess は「教育・研究・一次見積り用途で
  主要プロセス工程を一気通貫でシミュレートできる 3D プロセスシミュレータ」に到達した。
  MC 注入(化合物 BCA/チャネリング)、TED/クラスタ、OED、SPER、応力・SiGe 結合、GDS 入力、
  デバイス連携出力、決定的並列化と 100 本のテストは、この規模の実装としては例外的に整合している。
- 一方、**SProcess との差は「機能の有無」から「定量予測性と結合度」へ移った**。
  実務プロセス開発(レシピ転写・形状/Rs の定量予測)にそのまま代替使用できる段階ではなく、
  その距離を決めているのは (a) 校正データ、(b) 完全結合ソルバー+粘性流動酸化+自動リメッシュの
  三本柱、(c) フロントエンドの追随、である。
- 直近の一手としては**第 1 群(配線・更新)を先に片付ける**のが投資効率が最も高い。
  エンジンに既にある能力(AMG/BCSR/coarsen/oxidize_2d/silicide/…)を、変換器・デッキ・文書・
  主経路が使えていないことが、現時点の「見かけのギャップ」の相当部分を占めるからである。
