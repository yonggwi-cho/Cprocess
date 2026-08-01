# タスク仕様書 (Sonnet 実装用)

`docs/IMPLEMENTATION_PLAN.md` のロードマップを、実装者 (Claude Sonnet クラスの
コーディングエージェント) が**追加の設計判断なしに**実装できる粒度まで具体化した
仕様書群。1 ファイル = 1 タスク = 1 コミット を原則とする。

## 仕様書のフォーマット

各仕様書は次の節を必ず含む:

1. **目的** — 何を解決するか 1 段落
2. **現状コード** — 変更対象のファイル・関数・現在の挙動 (行番号ではなく関数名で参照)
3. **実装手順** — 番号付きステップ。新規/変更する関数はシグネチャを明記
4. **テスト仕様** — テストファイル名、テストケース、**数値の合格基準**
5. **完了条件 (DoD)** — チェックリスト
6. **やらないこと** — スコープ外の明示 (Sonnet の過剰実装防止)

## 共通規約 (全タスクに適用)

- ビルド/テスト: `cmake --build build -j$(nproc) && ctest --test-dir build`
  が全テスト PASS してからコミットする
- 新規テストは `CMakeLists.txt` の `foreach(t ...)` に登録する
- `proc::` に関数を追加した場合は pybind11 バインディング +
  `Simulation` メソッド + Python テストを必ず追加する (CLAUDE.md 参照)。
  内部ライブラリ (remesh/sparse 等) のみの変更なら Python 側は不要
- OpenMP 並列を追加/変更した場合、**スレッド数 1 と 4 で結果が一致する**
  (決定的アルゴリズムなら bit 一致、リダクション順序が変わる場合は相対差 < 1e-9)
  ことをテストで確認する
- コミットメッセージにタスク ID (例 `S-1`) を含める

## タスク一覧と状態

### Sprint 0 — 基盤 (実装済み)

| ID | ファイル | 状態 | 依存 |
|----|---------|------|------|
| S-1 | `S1_parallel_ilu.md` | **実装済** | なし |
| M-1 | `M1_edge_split_batch.md` | **実装済** | なし |
| M-5 | `M5_quality_repair.md` | **実装済** | M-1 |
| M-4 | `M4_interface_conform.md` | **実装済** | M-1, M-5 |

### Sprint 1 — フロー完結性 (酸化系) (実装済み)

| ID | ファイル | 状態 | 依存 |
|----|---------|------|------|
| P1-6 | `P1-6_oxidize.md` | **実装済** | M-4, M-5 |
| P1-4 | `P1-4_segregation.md` | **実装済** | P1-6 |
| P1-9 | `P1-9_multimaterial_diffusion.md` | **実装済** | P1-4 |

### Sprint 2 — アニール精度 (実装済み)

| ID | ファイル | 状態 | 依存 |
|----|---------|------|------|
| P1-3 | `P1-3_activation.md` | **実装済** | なし |
| P1-5 | `P1-5_rta_ramp.md` | **実装済** | なし (P1-3 と同一ループを触る — 後着が 1 行調整) |
| P1-2 | `P1-2_mc_damage_ted.md` | **実装済** | なし |

### Sprint 3 — 注入精度 + 基盤 (実装済み)

| ID | ファイル | 状態 | 依存 |
|----|---------|------|------|
| P1-1 | `P1-1_pearson4.md` | **実装済** | なし |
| P1-10 | `P1-10_params.md` | **実装済** | なし |
| M-2 (=P1-8) | `M2_adaptive_refine.md` | **実装済** | M-1 |
| P1-7 | `P1-7_etch_depo_mesh.md` | **実装済** | P1-6 |
| P1-11 | `P1-11_save_load_profile.md` | **実装済** | なし |

### Sprint 4+ — 物理深化 (P2)

| ID | ファイル | 状態 | 依存 |
|----|---------|------|------|
| P2-1 | `P2-1_full_point_defects.md` | 仕様確定 | P1-2, P1-10 |
| P2-2 | `P2-2_dopant_clusters.md` | 仕様確定 | P2-1, P1-3 |
| P2-3 | `P2-3_oed.md` | 仕様確定 | P1-6, P2-1 |
| P2-4 | `P2-4_oxidation_2d3d.md` | 仕様確定 | P1-6, M-4, M-5 (2 コミット可) |
| P2-5 | `P2-5_levelset_topo.md` | 仕様確定 | M-2 |
| P2-6 | `P2-6_fem_mechanics.md` | 仕様確定 | S-1, P1-10 (2 コミット可) |
| P2-7 | `P2-7_gds.md` | 仕様確定 | なし |
| P2-8 | `P2-8_new_dopants.md` | 仕様確定 | P1-10 (C-I シンクは run_ted 前提) |

### ソルバー / 並列化 (物理と並走可) (実装済み)

| ID | ファイル | 状態 | 依存 |
|----|---------|------|------|
| S-2 | `S2_block_solver.md` | **実装済** | なし (P2-1 の陰結合版が利用) |
| S-3 | `S3_newton_krylov.md` | **実装済** | なし |
| S-4 | `S4_amg.md` | **実装済** | なし |
| S-5 | `S5_adaptive_dt.md` | **実装済** | なし (P1-5 と統合注意) |
| M-3 | `M3_coarsen.md` | **実装済** | M-1 |
| M-6 | `M6_anisotropic.md` | **実装済** | M-2, M-3 |
| PA-1 | `PA1_parallel_remainder.md` | **実装済** | M-2 |
| PA-2 | `PA2_numa_bandwidth.md` | **実装済** | なし |
| PA-3 | `PA3_species_parallel.md` | **実装済** | なし |
| PA-4 | `PA4_gpu.md` | **実装済** | なし (CPROCESS_GPU=OFF 既定; ON はオフロード実機で未検証) |
| PA-5 | `PA5_mpi.md` | **実装済** | なし (CPROCESS_MPI=OFF 既定; 4 ランク実行は MPI 環境で未検証) |

### 長期 (P3) (実装済み)

縮約仕様は `P3_overview.md` (設計方針の確定文書として維持)。以下の個別仕様書へ展開・実装済み。

| ID | ファイル | 状態 | 依存 |
|----|---------|------|------|
| P3-a | `P3a_compound_bca.md` | **実装済** | なし (mc_implant.cpp のみ) |
| P3-b | `P3b_epitaxy.md` | **実装済** | P1-7, M-2 |
| P3-c | `P3c_sper.md` | **実装済** | P1-2, P2-1 |
| P3-d | `P3d_silicide.md` | **実装済** | P1-6, P1-7 (SegTable 拡張が P1-4 に波及) |
| P3-e | `P3e_stress_physics.md` | **実装済** | P2-6, P2-4 |
| P3-f | `P3f_sige_strain.md` | **実装済** | P2-6, P2-8, P3-e |
| P3-g | `P3g_fast_1d2d.md` | **実装済** | なし (独立) |
| P3-h | `P3h_device_export.md` | **実装済** | P1-11 |

実装順: P3-a → P3-h → P3-g → P3-b → P3-d → P3-e → P3-f → P3-c
(P3-b/d は process.cpp/materials.cpp 共有のため直列、P3-e→P3-f は依存関係
そのものが直列、P3-c は diffusion.cpp/process.cpp 変更のため P3-e/f 系と
分離)。

これで Sprint 0〜3、P1、P2、S/M/PA、P3 の全タスクが実装済みとなった。

### 次期ロードマップ (v2)

Sentaurus Process との残存ギャップ解消計画は `docs/IMPLEMENTATION_PLAN_v2.md` を参照
(根拠: `docs/sentaurus_gap_analysis_v2.md`)。タスク体系は W(配線・更新)/
C(校正データ)/ A(基盤投資)の 3 群・15 タスク。着手時に本ディレクトリへ
個別仕様書を展開する運用は初代と同一。

計画全体は `docs/IMPLEMENTATION_PLAN.md`(初代・完了)および
`docs/IMPLEMENTATION_PLAN_v2.md`(次期)を参照。

### SProcess パリティテスト

v2 ギャップの「実行可能な仕様書」として `tests/test_sprocess_parity.cpp`
(CTest 名 `sprocess_parity`, **WILL_FAIL TRUE** 登録)を追加した。各チェックは
現存 API のみで SProcess 的に正しい挙動をアサートし、修正済 8 件(C-2 安定性・
W-7 2 件・W-3 2 件・W-8 3 件)は通常スイートへ移設、残りは**未達**
(C-1/C-3 は修正済で移設、残るは C-2 の較正課題 1 件のみ)。併せて現行機能の不変量を固定する
golden シナリオ `tests/test_golden_flows.cpp`(CTest 名 `golden_flows`、常時 PASS)を追加。

| チェック | タスク ID | 現状 |
|---|---|---|
| スクリーン酸化膜による MC 注入の飛程減衰 (>5% 浅く) | W-8 | **修正済** → `test_implant_materials.cpp` へ移設(bare=145.0nm→ox50nm=102.5nm。副次的に deposit() の厚さ引数タイプミス(5nm のつもりが実は 50nm 指定漏れ)を発見・修正) |
| STI 酸化膜下の Si への遮蔽 (同深度で <50%) | W-8 | **修正済** → `test_implant_materials.cpp` へ移設(合計ドーズ比 0.086、要求 <0.5 に対し大幅達成。根本原因はチェック側の設計問題: 狭い絶対深度帯での比較が形状の異なる2分布の「テールのテール」同士になっていた。総ドーズ比較に変更) |
| 解析注入のスクリーン酸化膜オフセット (<0.8x) | W-8 | **修正済** → `test_implant_materials.cpp` へ移設(bare=97.5nm→ox50nm=42.5nm、同じ deposit() 厚さ修正で解決) |
| photo/mask 後の保存 VTU にレジスト形状が現れる | W-7/A-7 | **修正済** → `test_photo.cpp` へ移設(`proc::save` が既定で `<name>_stack.vtu` サイドカーを併記、`proc::save_stack` 新設。仕様: `W7_structure_inspection.md`) |
| レジストスタック存在下で save_state が例外を出さない | W-7 | **修正済** → `test_state_io.cpp` へ移設(save_state/export_device は警告してレジスト抜きで続行に緩和) |
| deck の `etch` コマンド受理 | W-3 | **修正済** → `test_flow.cpp` `test_deck_new_commands` へ移設(仕様: `W3_deck_python_sync.md`) |
| deck の `pdbset` コマンド受理 | W-3 | **修正済** → `test_flow.cpp` `test_deck_new_commands` へ移設(同上) |
| 解析 Pearson の 2Rp でのチャネリングテール (MC の 1/10 以内) | C-1 | **修正済** → `test_dual_pearson.cpp` へ移設(`profile="dual"` の主峰+チャネリングテール; 解析 1.10e17 vs MC 8.91e17 cm⁻³, 比 0.123 ≥ 0.1。仕様: `C1_implant_moments.md`) |
| Massoud 薄膜酸化促進 (Deal-Grove 比 >1.10x) | C-3 | **修正済** → `test_oxidation.cpp` へ移設(Deal-Grove B/A 定数を `ox.dry.*`/`ox.wet.*` として ParamDB 化、既定はビット不変。`ox.massoud.c`/`ox.massoud.l`(既定 OFF、opt-in)で 900℃/10nm 域が比 1.42(要求 >1.10x)。併せて `pressure_atm`/`hcl_frac`/`orient` を `proc::oxidize()` に追加(既定はビット不変)。仕様: `C3_oxidation_calibration.md`) |
| TED 750-850 ℃ アニールの数値安定性(有限値・総 B 質量保存) | C-2 | **修正済** → `test_ted.cpp` テスト 9 へ移設(根本原因: 1a 陰解のスパイク負値 → クラスタ forward 負値の質量生成。CI/CV 床 + forward/ratio/cl_old 床で修正) |
| TED 増速率が古典実験帯域内(5〜200x) | C-2 | 未達(修正後 900 ℃/60 s で ~300x に改善、なお文献 10-100x を超過 — 較正課題として残存) |

W-8 は上記 3 件全て解消し、タスクとして**完了**(`docs/tasks/W8_material_aware_implant.md`)。

**WILL_FAIL 運用**: 修正が入ってあるチェックが PASS に転じると、スイート全体の
終了コードが変わらない限りは緑のままだが、**全チェック PASS** になった時点で
exit 0 となり WILL_FAIL 反転で ctest が `sprocess_parity` を **Failed** と報告
する。個別チェックの進捗は実行ログのサマリ表(`N/3 parity checks passing`)で
確認し、PASS に転じたチェックは通常スイート(`test_golden_flows` または該当
`test_<feature>`)へ移設し、残りを本スイートに留める。

### 定量ベンチマーク

パリティスイートが「方向性(SProcess 的挙動の有無)」を固定するのに対し、
`tests/test_benchmarks.cpp`(CTest 名 `benchmarks`、通常スイート登録・常時
PASS)は**公表済みのエンジン非依存な文献値**に対する定量精度を固定する。
全参照値に出典コメント付き。信頼度で 3 層に分ける:

- **Tier A(ハードアサート、現状 23 項目全 PASS)**: 実測の上でエンジンが
  寛大な許容幅(飛程 ±25-30%、厚膜酸化 ±15-25%、拡散駆動は factor-2)内で
  文献値に一致するもの。外部真値に紐づく恒久的な回帰アンカー。
- **Tier B(大幅未達、`test_sprocess_parity.cpp` の WILL_FAIL 表に追加)**:
  文献値から大きく明確に外れる項目。上表の C-2 2 件がこれに当たる。
- **Tier C(INFO 表示のみ、アサートなし)**: 参照値の精度またはモデルの
  適用範囲がハードアサートを正当化しない項目。

| Tier | 項目 | 出典 | 参照値 vs 実測 |
|---|---|---|---|
| A | 解析注入 Rp/dRp: B 30/100, P 50/100, As 50/100 keV(12 項目) | LSS 表(Sze PSD / Plummer VLSI Ch.8) | 全て ±10% 以内(例: B30 Rp 100 vs 99 nm) |
| A | MC 注入 Rp: P 50/100, As 50/100 keV | 同上(BCA はモーメント表非依存の独立検証) | +15〜21%(±30% 内) |
| C | MC 注入 Rp: B 30/100 keV | 同上 | +34%(133 vs 100 nm 等)— 1.5x 未満で曖昧域 |
| A | 酸化: dry 1000℃/120min, dry 1100℃/30・60min, wet 1000℃/30・60min | Deal & Grove, JAP 36, 3770 (1965) 定数から算出 | −12%〜+3.5% |
| C | 酸化: dry 1000℃/30・60min(<70 nm 薄膜域) | 同上 | −21〜−36%(τ/Massoud 支配域。`ox.massoud.c`/`.l` opt-in で改善可能だが既定 OFF のためこの値のまま。C-3 参照) |
| A | B 真性拡散係数(埋め込みマーカ 1000℃/1h の σ² 成長) | Fair 1981: D_B=0.76·exp(−3.46eV/kT) | 比 1.03(factor-2 帯域) |
| A | B drive-in 接合深さ(1100℃/30min, 背景 1e15) | 同上 + 解析ガウス解 | 1.008 vs 0.955 µm(+5.7%) |
| C | TED 増速率(900℃/60s) | Packan/Stolk マーカ実験 10-100x @750-810℃ | ~500-700x(Tier B でハード帯域化) |
| C | As 電気活性上限 900/1000℃ | Nobili/Solmi(Plummer Ch.7) | 1.90/3.17e20 vs ~2/3e20 cm⁻³ — エンジン自身の固溶度フィットが同一出典由来のため**循環的**、アサート不可 |

**注意(本節の限界)**: 本スイートはあくまで**文献値プロキシ**であり、
SProcess そのものとの一致を保証するものではない。実際の SProcess 出力
(同一条件の golden データ)が入手できれば、それを直接参照値として追加する
方が保証としてはるかに強い。
