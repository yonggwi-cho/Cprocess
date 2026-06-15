# Sentaurus Process コマンド網羅マップ ↔ Cprocess

> Sentaurus Process (sprocess) で定義されている全コマンドを列挙し、Cprocess の
> 対応実装と、入力ファイル変換器 (`cprocess.sprocess`) での扱いをまとめる。
>
> 凡例:
> - ✅ **対応** — Cprocess に等価機能あり、変換器が自動翻訳
> - 🟡 **部分** — 近似/簡略化して翻訳（注記参照）
> - ⚠️ **無視** — 安全に読み飛ばし（警告コメントを出力）
> - ❌ **非対応** — 物理機能が未実装（変換器は警告コメント化）
>
> 作成日: 2026-06-15 / 出典: 末尾参照。

---

## 1. 構造・メッシュ定義

| Sprocess | 機能 | Cprocess 対応 | 状態 | 備考 |
|---|---|---|---|---|
| `line x\|y\|z location= spacing= tag=` | 1D グリッド線定義 | `mesh box` の extent/分割に集約 | 🟡 | line の最小/最大から box の範囲・分割数を推定 |
| `region <mat> xlo= xhi= ...` | 領域・材料定義 | `region material=` + box | 🟡 | Sprocess x(深さ)→Cprocess z |
| `init concentration= field=` | 基板初期化 | `init species= conc=` | ✅ | field=元素名→記号に変換 |
| `init tdr=<file>` | 構造ロード | `mesh gmsh file=` | ❌ | TDR 未対応。gmsh へ手動変換が必要 |
| `refinebox ... xrefine= yrefine=` | 適応細分化 | — | ⚠️ | 静的メッシュのため無視（box 分割で代替） |
| `grid remesh` / `grid set.*` | リメッシュ設定 | — | ⚠️ | 移動境界なしのため無視 |
| `transform reflect\|rotate\|cut` | 構造変換 | — | ❌ | 未対応 |
| `mater` / `material` | 材料パラメータ設定 | `region material=` | 🟡 | 既知材料(si/oxide/nitride/poly/gas)のみ |

## 2. プロセスステップ

| Sprocess | 機能 | Cprocess 対応 | 状態 | 備考 |
|---|---|---|---|---|
| `implant <species> dose= energy= tilt= rotation=` | イオン注入(解析) | `implant species= dose= energy= tilt= rotation=` | ✅ | 既定はガウシアン。元素名→記号変換 |
| `implant ... sentaurus.mc \| crystaltrim` | MC 注入 | `implant ... method=mc` | ✅ | チャネリング・ダメージ対応 |
| `diffuse temperature= time=` | アニール/拡散 | `diffuse temp= time=` | ✅ | ChargedFermi 相当 |
| `diffuse ... O2\|H2O\|gas_flow=` | 熱酸化 | — | ❌ | Deal-Grove 未実装。雰囲気は警告して拡散のみ実行 |
| `deposit material= thickness=\|rate= time= type=` | デポジション | — | ❌ | 移動境界未実装。警告コメント化 |
| `etch material= thickness=\|rate= time= type=` | エッチング | — | ❌ | 同上 |
| `strip <material>` | 層除去 | — | ⚠️ | 構造変化なしのため無視 |
| `mask name= left= right=` | マスク定義 | implant の `x1=x2=` 窓に転用 | 🟡 | 後続 implant の窓へ反映可 |
| `photo mask= thickness=` | レジスト塗布 | — | ⚠️ | 無視（mask 窓のみ保持） |
| `temp_ramp name= ...` | 温度傾斜 | `diffuse` を分割 | ❌ | 一定温度近似（平均温度で警告） |
| `anneal` | アニール(diffuse 別名) | `diffuse` | ✅ | diffuse と同義に翻訳 |
| `silicide` / `stress` | シリサイド/応力 | — | ❌ | 未対応 |
| `contact name= ...` | 電極定義 | — | ⚠️ | デバイス解析用、プロセスでは無視 |

## 3. パラメータ・制御

| Sprocess | 機能 | Cprocess 対応 | 状態 | 備考 |
|---|---|---|---|---|
| `pdbSet <Reg> <Sp> <Param> <val>` | パラメータDB上書き | — | ⚠️ | 設計方針により非再現。警告コメント化 |
| `pdbGet` / `pdbDelayDouble` | パラメータ取得 | — | ⚠️ | 同上 |
| `math numThreads= numThreadsMC=` | 並列数 | implant の `threads=` | 🟡 | MC スレッド数のみ転用 |
| `set <var> <val>` (Tcl) | 変数定義 | 変換器内で展開 | 🟡 | `$var` を字句置換 |
| `if/else/foreach` (Tcl) | 制御構文 | — | ❌ | 変換器は非対応（警告） |
| `fset` / `define` | マクロ | 字句置換 | 🟡 | 単純定義のみ |
| `exit` / `quit` | 終了 | `stop` | ✅ | |

## 4. 出力・後処理

| Sprocess | 機能 | Cprocess 対応 | 状態 | 備考 |
|---|---|---|---|---|
| `struct tdr=<file>` | 構造保存 | `save file=<file>.vtu` | 🟡 | TDR→VTU に拡張子変換 |
| `struct smesh=\|gmsh=` | メッシュ出力 | `save` | 🟡 | VTU で代替 |
| `select z=<expr>` | 後処理量選択 | — | ⚠️ | Python 側で `st.get_field` を使う想定 |
| `WritePlx <file>` | 1D プロファイル出力 | — | ⚠️ | Python で抽出（例: examples 参照） |
| `SetPlxList {...}` | 出力変数指定 | — | ⚠️ | 同上 |
| `layers` | 層厚レポート | — | ⚠️ | 構造変化なしのため無視 |
| `plot.1d` / `plot.2d` | 作図 | Plotly (Python) | 🟡 | `examples/visualize_3d.py` 参照 |
| `print` / `tclsel` | 値出力 | `print` | 🟡 | print に縮約 |

## 5. 種名・単位の対応

**元素名 → Cprocess 記号**（変換器が自動）:
`Boron→B, Phosphorus→P, Arsenic→As, Antimony→Sb`。
（As/B/P/Sb 以外は警告）。

**単位**: Sprocess は `<unit>` 角括弧で明示（例 `energy=20<keV>`、`dose=4e13<cm-2>`、
`time=30<min>`、`temperature=1000<C>`、`location=0.5<um>`）。
変換器は角括弧を剥がし Cprocess 接尾辞へ変換:
`<um>→um, <nm>→nm, <keV>→keV, <min>→min, <s>→s, <C>→C`。
単位省略時の Sprocess 既定（temperature=℃, length=µm, energy=keV, diffuse time=min,
dose=cm⁻²）を補う。

**座標系**: Sprocess は x をウェハ深さ方向（下向き正）とする。Cprocess は z を
深さとするため、変換器は **Sprocess x → Cprocess z** に対応付ける。

---

## 6. 変換器の使い方

`python/cprocess/sprocess.py` を参照。

```python
import cprocess as cp
from cprocess.sprocess import translate_file, run_sprocess

# (1) Sprocess .cmd を Cprocess デッキ文字列へ変換
deck = translate_file("nmos.cmd")
print(deck)            # 翻訳結果（未対応行は # [unsupported] コメント）

# (2) そのまま実行
st = run_sprocess("nmos.cmd")    # 内部で cp.run_deck を呼ぶ
B = st.get_field("B")
```

CLI:
```
python -m cprocess.sprocess nmos.cmd            # デッキを標準出力
python -m cprocess.sprocess nmos.cmd -o out.deck
python -m cprocess.sprocess nmos.cmd --run      # 変換して実行
```

### 変換の限界（重要）
- 移動境界を伴う `deposit/etch/oxidation` は構造を変えないため、結果は
  Sprocess と一致しない（注入・拡散の単体検証用と割り切る）。
- Tcl 制御構文（`if/foreach/expr`）は非対応。事前に展開された .cmd を想定。
- 3D box への近似のため 2D/3D の複雑構造は再現できない。

---

## 参照
- Sentaurus Process チュートリアル 1〜3, 15
  <https://ghzphy.github.io/Sentaurus_Training/sp/sp_1.html>,
  <https://ghzphy.github.io/Sentaurus_Training/sp/sp_2.html>,
  <https://ghzphy.github.io/Sentaurus_Training/sp/sp_3.html>
- Sentaurus Process User Guide (N-2017.09 公開版)
- 関連: `docs/sentaurus_reference.md`, `docs/sentaurus_gap_analysis.md`
