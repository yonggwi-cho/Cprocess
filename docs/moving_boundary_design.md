# 移動境界・適応メッシュ 設計ドキュメント

> Cprocess に非構造メッシュの移動境界（堆積 / エッチング / 熱酸化）と適応メッシュを
> 導入するためのアーキテクチャ設計。
>
> 作成日: 2026-06-15 / 対象ブランチ: `claude/quirky-bohr-7yfa3b`

---

## 0. 設計決定（確定事項）

| 論点 | 決定 | 根拠 |
|---|---|---|
| 界面発展方式 | **ALE 節点移動 ＋ レベルセット 併用** | 酸化は ALE（体積膨張・界面追従）、エッチ/デポはレベルセット（トポロジー変化に堅牢） |
| メッシュ生成器 | **内製（局所リメッシュ）** | 外部依存なし・ライセンス自由。edge flip/split/collapse + 平滑化で品質維持 |
| 次元スコープ | **最初から 3D** | 四面体非構造を直接対象。2D は 3D の薄板特殊ケースとして扱う |
| 力学結合 | フェーズ分割（初期は運動学のみ、後に粘弾性） | Deal-Grove 膜厚から界面速度を与える簡易結合で着手 |

既存の FVM 拡散ソルバ（`src/diffusion.cpp`）は面ベースで**任意の非構造四面体に対応済み**。
材料界面は「region タグが異なるセル間の内部面」として表現する。したがって本設計の
対象は **ジオメトリ発展とメッシュ再構成のみ**。

---

## 1. アーキテクチャ全体像

```
            ┌─────────────────────────────────────────────┐
            │              Deck / Python API              │
            │  deposit / etch / strip / oxidize / mask    │
            └───────────────┬─────────────────────────────┘
                            │
        ┌───────────────────┴────────────────────┐
        │            MovingBoundaryDriver          │
        │  (界面速度の計算と時間積分の統括)        │
        └───┬───────────────┬──────────────────┬───┘
            │               │                  │
   ┌────────▼──────┐ ┌──────▼────────┐ ┌───────▼────────┐
   │  ALE Mover    │ │  LevelSet     │ │  Remesher       │
   │ (節点移動+    │ │  (φ移流+      │ │ (局所操作:flip/ │
   │  メッシュ平滑)│ │  再初期化)    │ │  split/collapse)│
   └────────┬──────┘ └──────┬────────┘ └───────┬────────┘
            │               │                  │
        ┌───▼───────────────▼──────────────────▼───┐
        │              Mesh (既存) + 拡張           │
        │  nodes / cells / faces / cell_region      │
        │  + node_marker(界面/内部) + φ(level set)  │
        └───────────────┬───────────────────────────┘
                        │
            ┌───────────▼────────────┐
            │  FieldTransfer (保存的補間) │
            └────────────────────────┘
```

---

## 2. データ構造

### 2.1 Mesh 拡張（`include/cprocess/mesh.hpp`）

既存 `Mesh` を破壊せず、移動境界用フィールドを追加:

```cpp
struct Mesh {
    // ... 既存 ...
    // --- 移動境界用（未使用時は空） ---
    std::vector<int>    node_region_count;  // 節点が接する region 数（界面節点判定）
    std::vector<char>   node_on_surface;    // 露出表面節点フラグ
    std::vector<Vec3>   node_normal;        // 表面/界面の節点法線（外向き）
    std::vector<double> phi;                // レベルセット値（節点ベース、符号付き距離）
};
```

界面節点 = `node_region_count > 1` の節点。表面節点 = gas 領域に面する境界節点。

### 2.2 トポロジー補助構造（`include/cprocess/topology.hpp`, NEW）

局所リメッシュには「edge → 隣接 tet」「face → 隣接 tet」の高速参照が必須:

```cpp
struct MeshTopology {
    // edge (sorted node pair) -> incident cell ids
    std::unordered_map<uint64_t, std::vector<int>> edge_cells;
    // sorted node triple -> up to 2 incident cells（既存 face_lookup を活用）
    void build(const Mesh& m);
    uint64_t edge_key(int a, int b) const;  // a<b を 32bit ずつパック
};
```

### 2.3 材料・層情報（`include/cprocess/materials.hpp` 拡張）

```cpp
enum class Material { silicon, oxide, nitride, poly, photoresist, gas };
struct MaterialProps {
    double oxidant_diffusivity;  // 酸化剤拡散（Deal-Grove）
    double molar_volume_ratio;   // Si→SiO2 体積膨張比（2.2）
    // 後フェーズ: 粘度 η, ヤング率 E など
};
```

---

## 3. モジュール詳細

### 3.1 ALE Mover（`src/ale_mover.cpp`, NEW）

**役割**: 界面節点を物理速度で移動し、内部節点を平滑化して品質を保つ。

#### (a) 界面節点速度
- **酸化**: Deal-Grove の `dx/dt = B/(A+2x)` から界面法線速度 `v_n` を算出。
  Si 側は消費（後退）、SiO₂ 側は膨張（前進）。体積比 2.2 で節点変位を配分。
- **堆積**: 露出表面を法線方向に `rate·dt` 前進（等方）／指定方向（異方）。
- **エッチング**: 露出表面を法線方向に後退。マスク領域は固定。

法線は隣接境界面の面積加重平均:
```
n_node = normalize( Σ_f (face∋node) S_f )
```

#### (b) 内部節点の平滑化（メッシュ運動 PDE）
境界変位を内部へ滑らかに伝播。ラプラス平滑化（反復）:
```
x_i^{new} = (1-ω) x_i + ω · (Σ_{j∈neigh(i)} w_ij x_j) / Σ w_ij
w_ij = 1/|x_i - x_j|   (距離逆数重み)
```
境界節点は固定（Dirichlet）。要素反転チェック（符号付き体積 > 0）を毎反復で実施。
反転が出たら変位をサブステップ分割。

> 後フェーズ: 弾性アナロジー（線形弾性方程式を解いて変位場）に置換し、
> 大変形での品質を改善。さらに酸化は粘弾性流動へ。

### 3.2 LevelSet（`src/level_set.cpp`, NEW）

**役割**: エッチ/デポのトポロジー変化（アンダーカット・空隙・界面合体）を扱う。

#### (a) 表現
節点ベースの符号付き距離関数 `φ(x)`。`φ<0`=材料内, `φ>0`=外（gas）, `φ=0`=界面。

#### (b) 移流方程式
```
∂φ/∂t + v_n |∇φ| = 0     (法線方向速度 v_n で界面移動)
```
- エッチ: `v_n < 0`（材料後退）、デポ: `v_n > 0`。マスク下は `v_n=0`。
- 空間離散化: 節点勾配は最小二乗（既存 `gradients` を流用）。
  風上化（upwind）で安定化。
- 時間: 陽的、CFL `dt ≤ C·h_min/|v_n|`。

#### (c) 再初期化
移流で歪んだ φ を符号付き距離へ復元（Fast Marching / 反復 PDE）:
```
∂φ/∂τ + sign(φ₀)(|∇φ| − 1) = 0  → 定常まで
```

#### (d) 材料再構成 → メッシュへの反映
φ の符号からセルの region を更新。界面（φ=0）が tet を横切る場合、
**Remesher で界面に conforming な節点を挿入**（marching-tets 的な分割）し、
界面が面に乗るようにする。これが ALE との接続点。

### 3.3 Remesher（`src/remesh.cpp`, NEW）

**役割**: 内製の局所リメッシュ。品質劣化・界面挿入・粗密調整に対応。

#### 局所操作（3D 四面体）
| 操作 | 内容 | 用途 |
|---|---|---|
| **edge split** | エッジ中点に節点挿入、接する tet を分割 | 細密化、界面挿入 |
| **edge collapse** | エッジを潰し 2 節点を 1 に | 粗化、スリバー除去 |
| **2-3 / 3-2 flip** | 隣接 tet の対角面を入れ替え | 品質改善 |
| **4-4 flip** | 4 tet 構成の再配置 | スリバー対策 |
| **smoothing** | 節点位置最適化（Laplacian / 品質最適化） | 全体品質 |

#### 品質指標
```
q = 正規化体積/面積比 (radius ratio)  ∈ (0,1], 1=正四面体
スリバー = 体積≈0 だが辺長は有限の劣悪 tet → 最優先で除去
```
閾値 `q_min`（例 0.1）を下回る tet を検出し、flip/collapse/smoothing を適用。

#### 界面・物質境界の保存
- region 境界・外部境界の節点・エッジは collapse/flip で**潰さない**制約。
- 界面に conforming（界面が必ず face に乗る）を維持。

#### 実装順序
1. `MeshTopology` 構築（edge↔cell, face↔cell）
2. edge split（最も単純、界面挿入に必須）
3. Laplacian smoothing（既に ALE にあるものを共有）
4. edge collapse（境界保護付き）
5. 2-3/3-2 flip（品質改善）
6. スリバー除去（4-4 flip / collapse）

### 3.4 FieldTransfer（`src/field_transfer.cpp`, NEW）

**役割**: メッシュ変更時に濃度場を新メッシュへ転送。

- **ALE 節点移動のみ**（トポロジー不変）: 転送不要（同じ cell が動くだけ。
  ただし cell 体積が変わるので、**保存量 = 濃度×体積** を保つよう濃度を再スケール）。
- **局所リメッシュ / レベルセット再構成**（トポロジー変化）: 保存的補間。
  - 第1段階: 重心最近傍（単純・堅牢、軽度の数値拡散）
  - 第2段階: 供与体–受容体の体積重なり積分（質量厳密保存）

質量保存チェック: `Σ C_old·V_old ≈ Σ C_new·V_new`（テストで検証）。

### 3.5 MovingBoundaryDriver（`src/moving_boundary.cpp`, NEW）

統括ループ:
```
oxidize/deposit/etch(time, dt, ...):
  t = 0
  while t < time:
     1. 界面速度 v_n を物理から計算（Deal-Grove / rate）
     2. dt_sub = min(dt, CFL制限, 反転回避制限)
     3. ALE節点移動 or LevelSet移流 で界面を進める
     4. メッシュ品質チェック → 必要なら Remesher
     5. トポロジー変化があれば FieldTransfer
     6. cell_region / 幾何量を再計算（mesh.finalize 相当）
     t += dt_sub
```

---

## 4. デッキ / Python コマンド

```
deposit material=oxide thickness=0.02um [type=isotropic|anisotropic] [rate=]
etch    material=oxide thickness=0.01um [type=] [mask=m1]
strip   material=photoresist
mask    name=m1 x1=0um x2=0.2um y1=0um y2=0.5um   |   mask clear
photo   mask=m1 thickness=0.5um        # photoresist デポ + マスク設定
oxidize time=30min temp=1000C [wet=on|off] [dt=10s] [pressure=1]
```

既存コマンド（mesh/init/implant/diffuse/save/print）とは独立に追加。
移動境界系を未使用なら従来動作と完全互換。

Python へは pybind11 で `cp.deposit(...)` 等を公開、または `run_deck` 経由。

---

## 5. 物理モデル（数式は `docs/sentaurus_reference.md` 参照）

- **Deal-Grove 酸化**: `dx/dt = B/(A+2x)`、dry/wet の Arrhenius 係数、
  圧力・方位依存、Si消費＝酸化膜/2.2。界面で過剰格子間原子（OED, 後フェーズ）。
- **堆積/エッチ速度**: 等方（法線一様）／異方（指定方向）。マスクで遮蔽。

---

## 6. 段階的マイルストーン

| Phase | 内容 | 成果物 |
|---|---|---|
| **P0** | `MeshTopology` + edge split + Laplacian smoothing + 品質指標 | `topology.*`, `remesh.*`（split/smooth）, テスト |
| **P1** | ALE Mover（界面法線速度＋節点移動＋反転回避） + 平面デポ/エッチ | `ale_mover.*`, deck `deposit`/`etch`, テスト |
| **P2** | Deal-Grove 酸化（ALE、Si消費＋SiO₂膨張） | `oxidize` コマンド、解析解との比較テスト |
| **P3** | edge collapse + flip + スリバー除去（本格リメッシュ） | `remesh.*` 完成、品質テスト |
| **P4** | FieldTransfer 保存的補間（重なり積分） | `field_transfer.*`、質量保存テスト |
| **P5** | LevelSet（移流＋再初期化＋材料再構成） + マスク付きエッチ | `level_set.*`、トポロジー変化テスト |
| **P6** | 粘弾性力学結合（酸化の応力流動）、適応細密化の自動化 | 力学ソルバ、Sentaurus 相当 |

各 Phase は独立にビルド・テスト・コミット可能。P0〜P2 で「平面酸化・デポ・エッチ」が
動き、P3〜P5 で本格的な 3D 非構造移動境界、P6 で Sentaurus 同等を目指す。

---

## 7. リスクと対策

| リスク | 対策 |
|---|---|
| 3D 局所リメッシュのスリバー | radius-ratio 閾値検出 + 4-4 flip/collapse、最終手段で局所再生成 |
| 大変形での要素反転 | サブステップ分割、弾性アナロジー（P6）、反転検出で巻き戻し |
| 補間による数値拡散 | P4 で重なり積分の保存的補間に置換 |
| レベルセットと conforming メッシュの整合 | 界面挿入（marching-tets）で φ=0 を face に乗せる |
| 性能（毎ステップのリメッシュ） | 品質劣化時のみ局所操作、トポロジー不変時は転送スキップ |

---

## 8. 関連ドキュメント
- `docs/sentaurus_reference.md` — 物理モデル・数式
- `docs/sentaurus_gap_analysis.md` — Sentaurus との差分
- `docs/sprocess_command_map.md` — コマンド対応表
