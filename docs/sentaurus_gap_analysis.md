# Cprocess ↔ Synopsys Sentaurus Process 差分分析

> 業界標準である **Synopsys Sentaurus Process** に実装を寄せていくための
> ギャップ分析。物理モデル・数値計算手法・コマンド体系・データ構造の各観点で、
> Sentaurus の現行仕様と Cprocess の現状を対比する。
>
> 調査日: 2026-06-15 / 対象 Cprocess リビジョン: `claude/quirky-bohr-7yfa3b`

---

## 0. サマリ表

| 領域 | Sentaurus Process | Cprocess 現状 | ギャップ |
|---|---|---|---|
| 解析イオン注入 | Pearson-IV / dual-Pearson + Calibration Library | 単峰/双峰ガウシアン + LSS範囲表 | 中 |
| MC イオン注入 | Crystal-TRIM (BCA, 結晶/非晶質) | ZBL-BCA + Lindhard-Robinson チャネリング | 小 |
| 注入ダメージ | Kinchin-Pease + 蓄積非晶質化 + +1モデル | Kinchin-Pease + 非晶質化分率 | 中 |
| 拡散モデル | Fermi / ChargedFermi / Pair / ChargedPair / **five-stream(React)** | **ChargedFermi 相当**（中性/+/−/++ 成分） | 大 |
| 点欠陥 (I/V) | 陽に解く（注入・酸化・界面で生成消滅） | なし（平衡仮定） | 大 |
| 過渡増速拡散 TED | {311}/SMIC/クラスタと結合 | なし | 大 |
| ドーパント活性化 | None〜ComplexCluster の 9 モデル + 固溶度 | 固溶度クランプのみ（出力時） | 大 |
| 熱酸化 | Deal-Grove + 粘弾性応力 + 移動境界 | **なし** | 大 |
| エピ/エッチ/デポ | あり（構造編集 + レベルセット） | なし | 大 |
| 応力 | 粘弾性 / 塑性 | なし | 大 |
| メッシュ | 適応細分化, TDR | 静的非構造四面体 + box/gmsh | 中 |
| パラメータDB | `pdbSet`/`pdbGet` 階層DB | C++ ハードコード (`materials.cpp`) | 大 |
| コマンド言語 | Tcl 完全埋め込み | 独自 key=value パーサ | 中 |
| 入出力形式 | TDR (HDF5系) | VTU (ParaView) | 中 |
| 数値ソルバ | Newton 完全結合 + 直接/反復 | Picard + CG/BiCGSTAB/GMRES | 中 |
| 並列 | 共有/分散メモリ | OpenMP（MC・拡散） | 中 |

---

## 1. イオン注入

### 1.1 解析モデル
- **Sentaurus**: Pearson-IV 分布を基本とし、チャネリングテールを表現する
  dual-Pearson。Varian 等の実測を取り込んだ **Sentaurus Calibration Library**
  により高精度。マスク端の横方向散乱もモデル化。
- **Cprocess**: ガウシアン（`src/implant.cpp`）。範囲は LSS 風の粗い範囲表
  （`materials.cpp` の `range`、B/P/As/Sb のみ、対数補間）。`rp=`/`drp=` で上書き可。
  - **差分**: 分布形状が Pearson-IV でなくモーメント（歪度・尖度）を持たない。
    チャネリングテールを解析式では表現できない（MC では可）。

### 1.2 Monte Carlo (BCA)
- **Sentaurus**: Crystal-TRIM。結晶/非晶質ターゲット、電子的阻止能、
  ダメージ蓄積による動的非晶質化。
- **Cprocess**: `src/mc_implant.cpp`。ZBL 万有ポテンシャル + corteo 風 160×160
  対数散乱テーブル、Lindhard-Robinson 連続ストリングポテンシャルによる
  Si <100>/<110>/<111> チャネリング（13 軸）、Kinchin-Pease ダメージを
  OpenMP atomic で共有蓄積し `U_eff = U_max·(1−f_amor)` で動的非晶質化。
  - **差分（小）**: 物理は近接。Sentaurus のような実測校正テーブルは未整備。
    電子的阻止能は Lindhard-Scharff（E>500keV で警告）。

### 1.3 ダメージ／非晶質化
- **Sentaurus**: 蓄積ダメージ、"+1"/"effective-plus-factor" モデル、
  非晶質層の再結晶（SPER）。後段拡散の TED 初期条件に直結。
- **Cprocess**: 非晶質化分率 `f = damage_density / kNamorph`（`kNamorph=6.25e21`）。
  ただし**ダメージは後段拡散へ渡されない**（TED 未実装のため）。
  - **差分（中）**: ダメージ→欠陥→TED の連鎖が切れている。

---

## 2. 拡散

### 2.1 モデル階層
- **Sentaurus**: `pdbSet Silicon <species> DiffModel` で
  `Fermi → ChargedFermi → Pair → ChargedPair → React(five-stream)` を選択。
  five-stream は I, V, AI, AV, A の 5 連立を Newton で完全結合し、
  注入ダメージ起因の TED（過渡増速拡散）を再現。
- **Cprocess**: `src/diffusion.cpp` + `src/materials.cpp`。
  実効拡散係数を
  `D = D0·e^(−E0/kT) + D⁻·e^(−E⁻/kT)·(n/ni) + D⁼·(n/ni)² + D⁺·(n/ni)⁻¹`
  で表現（中性/単一負/二重負/単一正の荷電ペア）。
  これは Sentaurus の **ChargedFermi（平衡点欠陥）相当**。`n/ni` は全種の
  電荷中性から自己無撞着に決定（Picard 反復）。
  - **差分（大）**: 点欠陥 I/V を陽に解かない → TED・OED・エミッタ押し込み等の
    非平衡効果が出ない。{311}/クラスタによる過渡放出も無し。

### 2.2 物理係数（現状値）
`materials.cpp` の Fair 真空モデル定数（SUPREM 既定相当）:

| 種 | D0,E0(中性) | D⁻,E⁻ | D⁼,E⁼ | D⁺,E⁺ | 固溶度 pre,E |
|---|---|---|---|---|---|
| B | 0.037, 3.46 | – | – | 0.72, 3.46 | 9.25e22, 0.73 |
| P | 3.85, 3.66 | 4.44, 4.00 | – | – | 2.45e23, 0.62 |
| As | 0.066, 3.44 | 12.0, 4.05 | – | – | 1.3e23, 0.66 |
| Sb | 0.214, 3.65 | 15.0, 4.08 | – | – | 3.8e21, 0.56 |

`ni_si(T) = 3.87e16·T^1.5·e^(−0.605/kT)`。電界増速は `field_enh`（既定 on）。
  - **差分**: 係数は文献既定でハードコード。Sentaurus の `pdbSet` のような
    実行時上書き機構がない。

### 2.3 活性化・固溶度
- **Sentaurus**: `ActiveModel` に
  `None/Solid/Transient/Precipitation/Cluster/BIC/ChargedCluster/ComplexCluster/Equilibrium`。
  BIC（ボロン格子間クラスタ）等で非活性ドーパントを陽に追跡。
- **Cprocess**: 拡散中は全濃度を可動と仮定。活性濃度は**出力時のみ**
  固溶度でクランプ（`deck.cpp` `cmd_save`、`<species>_active`/`NetDoping`）。
  - **差分（大）**: クラスタ・析出の動力学なし。電気的活性は事後処理。

---

## 3. 熱酸化・構造変化

- **Sentaurus**: Deal-Grove（線形/放物速度定数）+ 粘弾性応力、酸化由来の
  界面欠陥生成（OED）、移動境界による体積膨張、全拡散ステップで native oxide
  1.5 nm を自動考慮。エピ成長・エッチング（レベルセット）・デポジション・
  CMP・構造編集。
- **Cprocess**: **未実装**。`region material=oxide` でラベル付けはできるが
  酸化成長・移動境界・応力は無い。
  - **差分（大）**: プロセス統合（注入→酸化→拡散の連成）が不可。

---

## 4. メッシュ・数値計算

### 4.1 離散化
- **Sentaurus**: ボックス法/有限要素、適応細分化、TDR 入出力。
- **Cprocess**: 非構造四面体 FVM（`src/mesh.cpp`）。box（`box_mesh.cpp`）/
  gmsh（`gmsh_reader.cpp`）読み込み。面の非直交補正（deferred correction、
  最小二乗勾配）、直交性メトリクス出力。
  - **差分（中）**: 静的メッシュ（適応細分化なし）。酸化の移動境界に未対応。

### 4.2 時間積分・非線形
- **Sentaurus**: 後退オイラー/TR-BDF、**Newton 法**で全種・全欠陥を完全結合。
- **Cprocess**: 後退オイラー（陰）+ **Picard（不動点）反復**で種間結合
  （`diffusion.cpp`、`max_picard`/`picard_tol`）。時間刻みは `dt` 指定 or `time/50`。
  - **差分（中）**: Newton ではなく Picard のため、強結合（five-stream）には
    収束性が不足。ヤコビアン未構築。

### 4.3 線形ソルバ（`src/sparse.cpp`, CSR）
- **実装済み**:
  - Jacobi 前処理 **CG**（対称系既定）
  - **BiCGSTAB**（CG 失敗時フォールバック、非対称対応）
  - **GMRES(m)**：Arnoldi + Givens 回転 QR + 右 Jacobi 前処理
    （Saad & Schultz 1986、`restart` 既定 30）
- **Sentaurus**: 直接法（疎 LU, PARDISO 系）と反復法、領域分割。
  - **差分（中）**: 直接法・ILU 等の高度前処理なし。Newton 結合系の
    ブロック構造に未対応。

### 4.4 並列
- **Cprocess**: OpenMP。MC はイオン並列（非晶質モードはビット再現、
  結晶モードは蓄積順依存で非決定）。拡散はセルループ並列。
- **Sentaurus**: 共有/分散メモリ。
  - **差分（中）**: MPI 分散なし。

---

## 5. ユーザインタフェース

### 5.1 コマンド言語
- **Sentaurus**: Tcl 完全埋め込み。`pdbSet/pdbGet`、`[Arrhenius pre Ea]`、
  変数・制御構文。例:
  ```tcl
  implant Boron dose=4e13 energy=20 tilt=7 rotation=0
  diffuse  temperature=1000 time=30
  diffuse  temperature=1200 time=120 H2O      ;# 湿式酸化
  pdbSet Silicon Boron Dstar {[expr 4.0*[pdbDelayDouble Silicon Boron Dstar]]}
  ```
- **Cprocess**: 独自 1 行 key=value パーサ（`src/deck.cpp`）。単位接尾辞
  （um/nm/min/keV/C 等）対応。例:
  ```
  mesh box xmax=0.5um ymax=0.5um zmax=1.2um nx=10 ny=10 nz=80
  init    species=P conc=1e15
  implant species=B energy=50keV dose=1e13 method=mc ions=80000 seed=7
  diffuse time=30min temp=1000C dt=30s
  ```
  対応コマンド: `mesh/region/init/implant/bc/diffuse/save/print/stop`。
  - **差分（中）**: 変数・式・制御構文・`pdbSet` 無し。種名は記号
    （B/P/As/Sb）で Sentaurus の元素名（Boron 等）と非互換。引数名も相違
    （`temp` vs `temperature`、`time=30min` vs `time=30`）。

### 5.2 パラメータデータベース
- **Sentaurus**: 材料×種×パラメータの階層 PDB を実行時に上書き可能。
- **Cprocess**: `materials.cpp` にハードコード（`kDopants`）。再コンパイル必須。
  - **差分（大）**: 校正フローに致命的。`pdbSet` 相当が必要。

### 5.3 Python / 入出力
- **Cprocess**: pybind11 で全 API を Python 公開（`python/`）。VTU 出力
  （ParaView）、Plotly 3D ブラウザ可視化。これは **Sentaurus に無い独自の強み**。
- **Sentaurus**: TDR（HDF5系）、Sentaurus Visual、SWB（ワークベンチ）。
  - **差分**: 形式互換なし（TDR 読み書き未対応）。

---

## 6. 推奨ロードマップ（段階的に標準へ寄せる）

優先度は「校正可能性 × 実装コスト × 物理的重要度」で評価。

1. **コマンド文法の Sentaurus 互換層**（低コスト・高効果）
   - 元素名エイリアス（Boron↔B）、`temperature=`/`time=`（無単位は ℃/秒）、
     `H2O`/`O2` フラグ受理（酸化フックの足場）。既存デッキとの後方互換維持。
2. **`pdbSet`/`pdbGet` パラメータDB**（中コスト・高効果）
   - `materials.cpp` の係数を実行時上書き可能な階層 KV ストアへ外部化。
     校正フローの前提。
3. **点欠陥 I/V + Pair 拡散**（高コスト・最重要物理）
   - I/V を陽に解き、Newton 結合へ移行。TED の土台。
4. **活性化・クラスタ（BIC/固溶析出）**（高コスト）
   - 拡散と結合した活性濃度の動力学。
5. **Deal-Grove 酸化 + 移動境界**（高コスト）
   - プロセス連成の完成。メッシュの動的更新が前提。

### 数値面の前提整備
- Picard → **Newton**（ヤコビアン構築、ブロック CSR）への移行は 3 以降で必須。
- 線形ソルバは ILU(0) 前処理 / 直接法フォールバックの追加を検討。
- メッシュ適応細分化は 5（酸化移動境界）で必要。

---

## 付録: 参照

- Synopsys TCAD / Process Simulation 製品ページ
  <https://www.synopsys.com/manufacturing/tcad/process-simulation.html>
- Sentaurus Process チュートリアル（モデル指定 `pdbSet`/`diffuse`/`implant`）
  <https://ghzphy.github.io/Sentaurus_Training/sp/sp_4.html>
- R.B. Fair, "Impurity Doping Processes in Silicon" (1981) — 拡散係数既定値
- Saad & Schultz (1986) — GMRES
