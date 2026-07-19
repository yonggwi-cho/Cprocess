# Cprocess 実装計画 v2 — Sentaurus Process 残存ギャップ解消ロードマップ

作成日: 2026-07-14 / 前提: HEAD `e3ef2bb`(全初代ロードマップ完了 + 総合評価 v2)
根拠: `docs/sentaurus_gap_analysis_v2.md` の残存ギャップ 14 項目(3 群)

本書は初代 `IMPLEMENTATION_PLAN.md` の後継である。運用は初代と同一:
本書で計画を確定し、**着手時に各タスクを `docs/tasks/` の個別仕様書
(6 節形式: 目的/現状コード/実装手順/テスト仕様/DoD/やらないこと)へ展開**してから
実装エージェントに渡す。共通規約(3 層アーキテクチャ、proc:: 関数への
pybind + Simulation メソッド + Python テスト義務、全テスト PASS でコミット、
仕様と実測が食い違う場合は計測に基づく逸脱をコード内に文書化)も初代を引き継ぐ。

タスク ID 規約: **W**(Wiring: 配線・更新)/ **C**(Calibration: 校正データ)/
**A**(Architecture: 基盤投資)。

---

## 0. 全体像と優先方針

| 群 | 性質 | 工数 | 効果 | スプリント |
|---|---|---|---|---|
| 第 1 群 W-1〜W-6 | 実装済み能力の配線・文書/変換器の追随 | 小 | 大(見かけのギャップの相当部分を解消) | Sprint 1〜2 |
| 第 2 群 C-1〜C-3 | 校正データ整備(定量予測性の主制約) | 中 | 大(日常使用の精度に直結) | Sprint 3 |
| 第 3 群 A-1〜A-6 | SProcess の三本柱への基盤投資 | 大 | 質的転換 | Sprint 4〜6 |

原則: **W 群を先に完了させる**。エンジンに既にある能力(AMG/BCSR/coarsen/
oxidize_2d/silicide/…)を変換器・デッキ・主経路が使えていないことが、
現時点の実効ギャップの最大成分だからである(v2 分析 §4 第 1 群)。

---

## Sprint 1 — 変換器・文書・API 同期(全て並列可)

### W-1: sprocess.py 変換器の全面更新【最優先】
- **現状**: `python/cprocess/sprocess.py` が 2026-06-15 時点のエンジン能力を前提。
  `diffuse+O2/H2O`→不活性アニールに劣化、`deposit/etch/silicide/temp_ramp`→unsupported、
  `photo/strip`→ignore、`pdbSet`→警告のみ、対応種 B/P/As/Sb のみ。変換器テスト無し。
- **やること**: (1) `diffuse` の雰囲気引数を `oxidize`/`oxidize_2d` へ変換
  (dry/wet 判定、時間・温度引き継ぎ)。(2) `deposit/etch`(rate×time は厚さへ換算)、
  `silicide`、`photo/strip`、`temp_ramp`(Python `ramp=` 経由)への対訳追加。
  (3) `pdbSet` → `set_param` のキー対訳表(対応可能なキーのみ、残りは警告)。
  (4) 8 種(+BF2→B への質量換算警告付き分解)対応。
  (5) `python/test_sprocess_converter.py` 新設(代表レシピ 5 本の変換結果ゴールデンテスト)。
- **DoD**: 変換器がエンジン実装済みコマンドを 1 つも unsupported にしない。
  代表 SProcess レシピ(LOCOS + well + S/D)が警告最小で実行可能なデッキ/Python に変換される。
- 依存: なし。規模: 中(1 タスク)。

### W-2: 陳腐化文書の更新
- **現状**: `sprocess_command_map.md` と初版 `sentaurus_gap_analysis.md` が
  実装済み機能(酸化/deposit/etch/silicide/stress/TED/PDB)を「未実装」と記載。
- **やること**: command_map を現行デッキ 16 コマンド + Python 37 メソッドに合わせて全面改訂
  (W-1 の対訳表と整合させる)。初版 gap_analysis の冒頭に「v2 (`sentaurus_gap_analysis_v2.md`)
  で置換済み」の注記を付し、サマリ表を現状に更新。README のユーザ向け導線
  (どの文書が最新か)を明記。
- **DoD**: 全文書の機能記載が HEAD の実装と一致。git grep で「未実装」と実装の矛盾ゼロ。
- 依存: W-1 と同時進行可(対訳表のみ W-1 完了後に反映)。規模: 小。

### W-3: デッキと Python API の機能同期
- **現状**: デッキ 16 コマンド vs Python 37 メソッド。etch/oxidize_2d/sper/mechanics/
  refine/mask_polygon/ramp/set_param/save_state/export_device がデッキから呼べない。
  `Simulation.bbox()` が µm 系 API の中で唯一 cm を返す(simulation.py:568)。
- **やること**: (1) deck.cpp に `etch`(blanket/material/polygon)、`oxidize2d`、`sper`、
  `mechanics`、`refine`、`ramp`(diffuse のオプション)、`pdbset`(ParamDB 設定)、
  `save_state`/`load_state`、`export_device` を追加。単位接尾辞・行番号エラーの既存規約踏襲。
  (2) `bbox()` を µm 返しに修正(docstring と挙動の一致、既存呼出し箇所を確認して同時修正)。
  (3) ParamDB に未知キー警告(`log` へ warn、動作は従来どおり不活性)を追加。
- **DoD**: デッキで NMOS LOCOS フロー全体(注入〜シリサイド〜export)が書ける。
  examples/ にデッキ版フルフローを 1 本追加し test_flow 相当で検証。
- 依存: なし。規模: 中。

### W-7: 途中構造の検査機能【外部テスト指摘による追加】
- **現状**: レジストは本体メッシュと別の隠しスタック(`SimState::stack`)にあり、
  `save()` は本体メッシュのみ出力、`save_state()`/`export_device()` はレジストが
  あると throw。**photo/mask_polygon の結果形状を確認する手段が皆無**
  (根本分析: `docs/structure_model_root_cause.md`)。
- **やること**: (1) `proc::save_stack(st, path)` — stack メッシュを材料インデックス
  (Si/resist/open)付き VTU 出力(既存 write_vtu 流用)。(2) `save()` に
  `include_stack` オプション(レジスト存在時 `<name>_stack.vtu` を併記)。
  (3) Python `Simulation.resist_mask()`(セル中心 + 材料の numpy 配列)。
  (4) save_state/export_device の throw を「警告してレジスト抜きで続行」に
  緩和(実装時判断: パリティチェックが素の呼び出しの非 throw を要求し、
  photo/mask 状態は安価に再構築できるため、`force=` 待避なしで既定を
  警告+続行に変更。詳細: `docs/tasks/W7_structure_inspection.md`)。
- **DoD**: 複雑ポリゴンマスク(GDS 読込含む)を photo→mask_polygon→save_stack で
  ParaView 確認できる。pybind + Simulation + Python テスト三点セット。
- 依存: なし。規模: 小。

### W-8: メッシュ材料を考慮した注入輸送【外部テスト指摘による追加】
- **現状**: MC 輸送エンジンは多材料対応済み(P3-a、oxide/nitride 化合物 BCA まで)
  だが、一般経路は材料テーブルを渡さず**メッシュ全体を結晶 Si として輸送**。
  STI/スクリーン酸化膜/窒化膜マスク越しの注入が物理的に誤る。解析注入も
  上層材料の深さオフセットなし。多材料テーブルはレジストスタック経路のみ
  ({Si, resist, 真空} の 3 種限定)。
- **やること**: (1) `material_ids(st)`(既存)→ TargetMaterial テーブル
  ({Si, SiO2, Si3N4, poly≈Si, gas=真空})写像を構築し、一般 MC 経路で常時
  `cell_material` を渡す。(2) レジストスタック経路のテーブルにも oxide/nitride を
  追加(スタック下の実材料を反映)。(3) 解析注入にカラム毎スクリーニング補正
  (上層材料厚の実効 Si 換算)。(4) 全域 Si メッシュでは従来とビット一致を保証
  (single_material 高速パスの維持)。
- **DoD**: スクリーン酸化膜 20 nm 越し B 30 keV の Rp シフトが MC/解析で整合。
  STI 構造(酸化膜埋込トレンチ)への well 注入で酸化膜下の分布が bare Si と
  有意に異なることを検証。golden flow テスト(下記 GF)①③が PASS。
- 依存: なし(W-7 と並列可)。規模: 中。

### GF: golden flow シナリオテスト新設【再発防止・W-8 と同時】
- 実プロセスフロー横断の統合テスト 3 本を tests/test_golden_flows.cpp として新設:
  ① STI(トレンチエッチ+酸化膜埋込)→ well 注入 → RTA、
  ② LOCOS + poly ゲート + レジスト 2 回(S/D 注入)、
  ③ スクリーン酸化膜越し注入 + スパイク RTA + Rs 相当量の検証。
  以後の**全タスクの DoD に golden flow PASS を含める**(共通規約 7 として追記)。
- 依存: W-8(①③が W-8 の挙動を前提)。規模: 中。

### W-4: 既定スレッド数抑制の解除
- **現状**: `ensure_sane_thread_count()`(sparse.cpp)が OMP_NUM_THREADS 未設定時に
  既定 1 スレッドへ抑制、`run_ted` は `OmpThreadGuard(1)` で全体シングルスレッド化。
  CI 環境のスケジューラ問題対策だが実マルチコアの既定性能を殺している。
- **やること**: (1) 環境変数 `CPROCESS_THREADS`(または ParamDB `omp.threads`)で明示制御。
  (2) 既定は「物理コア数と 8 の小さい方」等の安全な複数スレッドに変更し、
  オーバーサブスクライブ検知(起動時の軽量ベンチ)で 1 に落ちる自動フォールバックを検討。
  (3) run_ted の OmpThreadGuard(1) を、当時の実測ハング原因(反応サブサイクルの
  小粒度並列)を再計測した上で、問題のあるループのみのガードに縮小。
  (4) 変更前後で全テストのビット一致/決定性テスト(topology_parallel 等)を必ず確認。
- **DoD**: 4 コア環境で ted/oed/ox2d 系テストの実測 wall time が既定設定で短縮。
  決定性テスト全 PASS。
- 依存: なし。ただし性能計測を伴うため専用タスクとして単独実行。規模: 中。

---

## Sprint 2 — 実装済み部品の主経路統合(W-5, W-6 は直列推奨: 双方 diffusion.cpp)

### W-5: AMG・ブロック解法の主経路統合
- **現状**: `cg_amg`(2 レベル SA-AMG)と `BCSR`+`bicgstab_bjacobi`(S-2/S-4)は
  テストのみが使用。`DiffusionSolver` は常に ILU0、`run_ted` は I/V/ドーパント逐次分離。
- **やること**: (1) `DiffuseOpts::precond = "ilu0"|"amg"|"auto"` を追加し、
  `solve_permuted()` で選択(auto: nc 閾値で切替、閾値は実測で決定)。
  (2) run_ted の反応サブサイクル内 I-V 結合(双分子再結合 + {311})を
  2×2 BCSR ブロック陰解に置換するオプション(`DiffuseOpts::coupled_iv`、既定 off)。
  既定 off で全既存テストのビット不変を保証、on で剛性の強いケース
  (高ドーズ TED)の許容 dt 拡大を実測で示す。
- **DoD**: AMG 選択時に 10⁶ セル級で ILU0 比の反復数半減(S-4 テスト基準の主経路版)。
  coupled_iv on で test_ted 相当ケースが同精度・より少ないサブサイクルで PASS。
- 依存: W-4(スレッド既定変更後に性能基準を取り直すため)。規模: 大。

### W-6: coarsen の配線 + 自動リメッシュの最小版
- **現状**: `coarsen`/`select_coarsen_edges`(M-3)は proc::/Python 未配線。
  refine は手動・操作単位。プロセス中の自動適応リメッシュなし。
- **やること**: (1) `proc::coarsen(st, species, rel_grad_thresh, max_fraction, log)` +
  pybind + `Simulation.coarsen()` + Python テスト(CLAUDE.md 規約)。
  (2) 最小の自動化: `DiffuseOpts::auto_remesh`(既定 off)で、diffuse の
  採択ステップ N 回ごとに refine(高勾配)→ coarsen(低勾配)→ repair_quality を
  1 パス実行するフック。質量保存(<1e-9)と min_q>0.05 を毎パス assert。
  (3) refinebox 相当の領域限定 refine(bbox 指定)を refine に追加。
- **DoD**: TED ケースで auto_remesh on がセル数を有界に保ちつつ固定細メッシュと
  ピーク/ドーズ相対差 <2%。既定 off で全既存テストのビット不変。
- 依存: W-5 完了後(diffusion.cpp 共有)。規模: 大。

---

## Sprint 3 — 校正データ整備(C-1〜C-3 並列可、いずれも materials.cpp 中心だが行区画が分離)

### C-1: 注入モーメント表の拡充 + dual-Pearson
- **現状**: 8 種×5–8 点(10–200 keV、教科書級 LSS 風)、範囲外クランプ、単峰 Pearson-IV。
  解析注入は日常使用の大半を占め、実効精度の最大の制約(v2 分析 物理 Top1)。
- **やること**: (1) SRIM/文献(Janson 等)に基づく B/P/As/Sb/In/Ge の
  1 keV〜3 MeV モーメント表(≥15 点/種、γ/β 込み)への差し替え。出典をコード内に明記。
  (2) dual-Pearson(主峰 + チャネリングテール峰、9 パラメータ)を
  `profile="dual"` として追加。テール峰パラメータは Si <100> 無傾角の代表条件で校正し、
  ParamDB キー(`<Sym>.dp.*`)で上書き可能に。(3) 検証: 代表条件(B 5/20/80 keV 等)で
  MC 結果とテール一致(ドーズ加重 L2 <15%)、SRIM Rp と <5%。
- **DoD**: 上記検証テスト + 既存 pearson/implant テストの期待値更新(出典コメント付き)。
- 依存: なし。規模: 大(データ整備が主)。

### C-2: クラスタ/活性化パラメータの文献校正 + P/Sb/In 拡張
- **現状**: B(BIC)/As(As4V)のみ、Eb はテスト合わせ校正(B 2.7 / As 2.6 eV と明記)。
  P/Sb/In クラスタなし。Transient/析出モデルなし。
- **やること**: (1) B/As の kf/ν0/Eb を文献(Pelaz, Solmi 等)ベースで再校正し、
  代表 RTA 条件(950–1050 ℃ スパイク)の活性化率が実測レンジに入ることを検証テスト化。
  (2) P に P-V/P4V 型実効クラスタ、Sb/In に固溶度クランプ強化(In の低固溶度)を追加。
  (3) SIMS/Rs 突合フローの整備: `Simulation.sheet_resistance()`(Irvin 曲線近似)と
  `junction_depth()` を追加(W 群の抽出機能とも整合)— pybind + Python テスト込み。
- **DoD**: B 1e15/20keV + 1000 ℃ RTA の Rs が文献実測 ±30%。既存 clusters/activation
  テストは校正値更新後も物理基準(質量保存・単調性)で PASS。
- 依存: なし。規模: 中〜大。

### C-3: 酸化の薄膜補正・雰囲気依存 + Deal-Grove 係数の ParamDB 化
- **現状**: DG 係数ハードコード、Massoud 補正なし、1 atm/HCl なし/<100> 固定。
- **やること**: (1) `ox.dry.b0/be/a0/ae`(wet 同様)キーで DG 係数を ParamDB 化
  (既定値は現行定数、ビット不変)。(2) Massoud 薄膜増速項
  (C·exp(−x/L)、<30 nm 域)を `ox.massoud.*` キー付きで追加(既定 on/off は実測で決定、
  off 時ビット不変)。(3) 圧力スケーリング(B∝P, B/A∝P^0.75)と HCl 増速係数、
  方位係数(<111>=1.68×<100>)を `oxidize(pressure=, hcl=, orient=)` として
  proc::/pybind/Simulation/デッキに追加。
- **DoD**: 薄膜域(10 nm)の成長が Massoud 文献曲線 ±10%。厚膜域は現行 DG とビット一致
  (補正off)。全酸化系テスト PASS。
- 依存: なし(W-3 のデッキ拡張と引数名を揃える)。規模: 中。

---

## Sprint 4〜6 — アーキテクチャ投資(A 群、順序に依存関係あり)

### A-1: 疎直接法の導入(Sprint 4)
- **内容**: 外部ライブラリ(Eigen SparseLU / SuiteSparse UMFPACK のいずれか、
  CMake オプションで opt-in)を `SolveResult solve_direct(CSR,...)` として統合し、
  反復法不収束時の自動フォールバックに配線(現行の即 throw を置換)。
  ライセンス・ビルド依存を仕様書で確定してから着手。
- **理由**: 悪条件系のロバスト性の根源。A-2(完全結合 Newton)の前提。
- 依存: なし。規模: 中。

### A-2: 完全結合多種 Newton(Sprint 4〜5、A-1 後)
- **内容**: run_ted の I/V/ドーパント(+クラスタ)を BCSR ブロックヤコビアン
  (解析微分: 拡散項は既存 assemble、反応項は質量作用の閉形式微分)で
  完全結合 Newton 化(`DiffuseOpts::fully_coupled`、既定 off)。
  線形解は bicgstab_bjacobi、フォールバックは A-1 の直接法。
  S-3 JFNK の「アフィン部分問題」限界(コード内文書化済み)の根本解。
- **理由**: 五流(React)級モデル(A-5)の必須基盤(reference §3.3「Picard では収束しない」)。
- 依存: A-1、W-5(BCSR 配線)。規模: 特大。2〜3 コミット分割(組立/求解/検証)。

### A-3: 2D/3D 酸化の粘性流動 + 一般移動境界(Sprint 5)
- **内容**: 酸化膜を Stokes 近似(非圧縮粘性流、粘性は温度依存 η(T))で解き、
  界面成長速度場から ALE ノード移動(既存 ale_mover 流用)+ 自動品質修復で
  一般メッシュの移動境界酸化を実現。現行カラム 1D `oxidize_2d` は高速パスとして残す。
  FEM(fem.cpp)の Stokes 拡張(混合要素 or 安定化 P1-P1)が主工数。
- **理由**: LOCOS/STI 形状の定量化(v2 分析 物理 Top3)。`ox2d.nitride_leak`
  チューニング係数からの脱却。
- 依存: A-1(Stokes 系は直接法が事実上必須)。規模: 特大。3 コミット分割
  (Stokes ソルバー/界面速度・ALE 結合/検証)。

### A-4: 自動適応リメッシュの本格化 + 界面適合(Sprint 5〜6)
- **内容**: W-6 の最小版を発展させ、(1) 誤差指標(勾配 + ヘッシアン近似)による
  refine/coarsen の自動オーケストレーション、(2) 界面近傍の能動的適合
  (界面ノードのスナップ + カットセル的 retag の廃止)、(3) 保存的フィールド転写
  (transfer_field_nearest の O(N²) 最近傍 → 交差体積ベース)への置換。
- 依存: W-6。規模: 特大。
- **注**: Kuhn-tet 非直交問題(v2 数値 Top4)はここで **Voronoi/ボックス法双対格子の
  導入判断**を行う。導入する場合は別タスク A-4b として仕様展開(diffusion.cpp の
  build/assemble の離散化置換、全テスト期待値の再基準化を伴う最大級の変更)。

### A-5: 五流物理(荷電点欠陥 + AI/AV ペア)(Sprint 6、A-2 後)
- **内容**: I/V の荷電状態(I⁻/I⁰/I⁺、V²⁻…)と Fermi 準位依存分配、
  AI/AV ペア場の陽解(5 種連立)、Frenkel 質量作用。A-2 の完全結合 Newton 上に実装。
  既存 3 場モデルは `ted.model = "effective"|"five_stream"` で選択制。
- 依存: A-2。規模: 特大。
- 同時に検討: 注入間累積非晶質化の BCA 反映(damage 場→次回 MC の初期結晶状態)、
  分子注入 BF2(質量分配 + 実効エネルギー)— いずれも独立の中規模タスクとして分離可。

### A-7: 構造モデル統一(structure-first 化)(Sprint 4 冒頭、A-3/A-4 より先)
- **内容**: 構造表現の分裂(本体メッシュ/レジスト隠しスタック/layer_stack/
  レベルセット一時場/カラム高さ — `docs/structure_model_root_cause.md` §1)を解消する。
  第 1 段: レジストを本体メッシュの材料(MatId に resist 追加)として統合し、
  photo/mask を deposit/etch の特殊形へ再定義(既存 API は互換ラッパ維持、
  MC は W-8 の多材料輸送でそのまま動く)。表現 2(隠しスタック)を廃止。
  第 2 段: レベルセット φ とカラム高さを「mesh+材料から導出されるキャッシュ」に
  格下げし、正本を常に単一化。
- **理由**: SProcess の中核アーキテクチャ不変量(全工程が単一構造を共有)への
  収斂。A-3(粘性流動酸化)・A-4(界面適合)を分裂した表現の上に建てると
  表現がさらに増殖するため、**A-3/A-4 より先に第 1 段を完了させる**。
- 依存: W-7, W-8(是正の前提)。規模: 大。2 コミット分割(第 1 段/第 2 段)。

### A-6: デバイス連携の実用化(Sprint 6、独立)
- **内容**: (1) `contact` コマンド(名前 + 面パッチ/領域指定)を SimState に追加し
  export_device の meta.json に出力。(2) TDR は非公開フォーマットのため直接対応は
  見送り、代替として SDevice 側で meta.json + VTU を読む変換スクリプト
  (リファレンス実装)を tools/ に整備。(3) 抽出強化: `layers()`(層厚レポート)、
  WritePlx 互換 1D 出力。C-2 の Xj/Rs と合わせ SProcess extract 相当を概ね充足。
- 依存: なし(C-2 の抽出と整合)。規模: 中。

---

## 実行順序まとめ(依存グラフ)

```
Sprint 1:  W-7 ∥ W-8 ∥ W-1 ∥ W-2 ∥ W-3 ∥ W-4   (W-7/W-8 最優先; W-2 の対訳表のみ W-1 後)
Sprint 1':  GF(golden flow テスト; W-8 後すぐ)
Sprint 2:  W-5 → W-6                        (直列: diffusion.cpp 共有)
Sprint 3:  C-1 ∥ C-2 ∥ C-3                 (並列可; W-3 と引数規約を整合)
Sprint 4:  A-7 第1段(レジスト統合) → A-1 → A-2(開始) ∥ A-6
Sprint 5:  A-2(完了) → A-3 ∥ A-4(W-6, A-7 後)
Sprint 6:  A-5(A-2 後) ∥ A-4 継続・A-4b 判断 ∥ A-7 第2段
```

## 共通規約(初代から継承 + 追加)

1. 着手時に `docs/tasks/` へ 6 節個別仕様書を展開してから実装(W 群は 1 仕様 1 タスク、
   A 群は複数コミット分割を仕様書内で明示)。
2. proc:: 新関数は pybind + Simulation + Python テストの三点セット必須(CLAUDE.md)。
   **Python テストは定義追加と `__main__` 呼出しリスト登録の両方を確認**
   (初代で 2 度発生した登録漏れの再発防止)。
3. 既定 off のオプション追加時は off 経路のビット不変をテストで保証(S-5/PA-3/P3-e 方式)。
4. 仕様の数値基準が実測と食い違う場合は、計測プログラムによる根拠と共にコード内へ
   逸脱を文書化(既存コミット 777531e〜7055e35 の確立パターン)。
5. 校正系(C 群)は係数の出典(文献/SRIM 条件)をコードコメントに必須記載。
6. フルテストスイート(現 49 本 + Python 51)は今後も増加するため、
   長時間テスト(ted/oed/ox2d/newton)の分離実行ラベル化を W-2 で合わせて導入検討。
7. **golden flow ゲート**(GF 完了後): 全タスクの DoD に golden flow 3 本の PASS を
   含める。工程横断の構造不整合をタスク単位検証の外に漏らさない
   (`docs/structure_model_root_cause.md` §3 Why-3 への恒久対策)。
8. **アーキテクチャ不変量**: 全工程は SimState の単一構造表現を読み書きする。
   工程私有の構造表現を新設する場合は ADR で例外理由と統合計画を必須記載
   (同レポート §4.2-1; CLAUDE.md にも追記済み)。
9. **配線負債レジスタ**: 「実装済みだが主経路未統合」のコンポーネントを
   docs/tasks/README.md の常設表で管理し、タスク完了時に新規負債を登録
   (現在の負債: AMG, BCSR, coarsen, GPU CG → W-5/W-6 で解消予定;
   多材料 MC 輸送 → W-8 で解消予定)。
10. 新工程機能の仕様書に「エンジニアはこの工程の結果をどう確認するか」節
   (可視化・抽出・ログ)を必須化(同レポート §4.2-5)。

## スコープ外(本計画では扱わない)

- SiC/Ge/III-V 基板対応(需要確定まで保留 — 材料表・注入表・酸化則の全面追加が必要)
- CMP、フォトリソの光学モデル(露光シミュレーション)
- TDR バイナリの直接読み書き(非公開フォーマット; A-6 の代替経路で対応)
- GPU 本格対応(cg_jacobi_gpu の拡散接続は W-5 の precond 選択肢設計時に判断のみ行う)
- MPI の拡散ソルバー分散(METIS 分割; PA-5 仕様書の将来項目のまま)
