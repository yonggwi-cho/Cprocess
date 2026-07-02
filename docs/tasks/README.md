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

### Sprint 1 — フロー完結性 (酸化系)

| ID | ファイル | 状態 | 依存 |
|----|---------|------|------|
| P1-6 | `P1-6_oxidize.md` | 仕様確定 | M-4, M-5 |
| P1-4 | `P1-4_segregation.md` | 仕様確定 | P1-6 |
| P1-9 | `P1-9_multimaterial_diffusion.md` | 仕様確定 | P1-4 |

### Sprint 2 — アニール精度

| ID | ファイル | 状態 | 依存 |
|----|---------|------|------|
| P1-3 | `P1-3_activation.md` | 仕様確定 | なし |
| P1-5 | `P1-5_rta_ramp.md` | 仕様確定 | なし (P1-3 と同一ループを触る — 後着が 1 行調整) |
| P1-2 | `P1-2_mc_damage_ted.md` | 仕様確定 | なし |

### Sprint 3 — 注入精度 + 基盤

| ID | ファイル | 状態 | 依存 |
|----|---------|------|------|
| P1-1 | `P1-1_pearson4.md` | 仕様確定 | なし |
| P1-10 | `P1-10_params.md` | 仕様確定 | なし |
| M-2 (=P1-8) | `M2_adaptive_refine.md` | 仕様確定 | M-1 |
| P1-7 | `P1-7_etch_depo_mesh.md` | 仕様確定 | P1-6 |
| P1-11 | `P1-11_save_load_profile.md` | 仕様確定 | なし |

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

### ソルバー / 並列化 (物理と並走可)

| ID | ファイル | 状態 | 依存 |
|----|---------|------|------|
| S-2 | `S2_block_solver.md` | 仕様確定 | なし (P2-1 の陰結合版が利用) |
| S-3 | `S3_newton_krylov.md` | 仕様確定 | なし |
| S-4 | `S4_amg.md` | 仕様確定 | なし |
| S-5 | `S5_adaptive_dt.md` | 仕様確定 | なし (P1-5 と統合注意) |
| M-3 | `M3_coarsen.md` | 仕様確定 | M-1 |
| M-6 | `M6_anisotropic.md` | 仕様確定 | M-2, M-3 |
| PA-1 | `PA1_parallel_remainder.md` | 仕様確定 | M-2 |
| PA-2 | `PA2_numa_bandwidth.md` | 仕様確定 | なし |
| PA-3 | `PA3_species_parallel.md` | 仕様確定 | なし |
| PA-4 | `PA4_gpu.md` | 仕様確定 | なし (長期) |
| PA-5 | `PA5_mpi.md` | 仕様確定 | なし (長期) |

### 長期 (P3)

| ID | ファイル | 状態 |
|----|---------|------|
| P3-a〜h | `P3_overview.md` | 縮約仕様確定 (着手時に個別展開) |

計画全体は `docs/IMPLEMENTATION_PLAN.md` を参照。
