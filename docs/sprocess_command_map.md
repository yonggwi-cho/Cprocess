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
> 作成日: 2026-06-15 / 更新: 2026-08-01 (W-1: 変換器全面更新に追随)
> 出典: 末尾参照。

---

## 1. 構造・メッシュ定義

| Sprocess | 機能 | Cprocess 対応 | 状態 | 備考 |
|---|---|---|---|---|
| `line x\|y\|z location= spacing= tag=` | 1D グリッド線定義 | `mesh box` の extent/分割に集約 | 🟡 | line の最小/最大から box の範囲・分割数を推定 |
| `region <mat> xlo= xhi= ...` | 領域・材料定義 | `region material=` + box | 🟡 | Sprocess x(深さ)→Cprocess z |
| `init concentration= field=` | 基板初期化 | `init species= conc=` | ✅ | field=元素名→記号に変換 |
| `init tdr=<file>` | 構造ロード | `mesh gmsh file=` | ❌ | TDR 未対応。gmsh へ手動変換が必要 |
| `refinebox ... xrefine= yrefine=` | 適応細分化 | `refine species= thresh= passes= axis=` | 🟡 | デッキに `refine` が存在するが `refinebox` の矩形指定とは意味が異なるため無視のまま。手動で `refine` 行を追加すること |
| `grid remesh` / `grid set.*` | リメッシュ設定 | — | ⚠️ | 移動境界なしのため無視 |
| `transform reflect\|rotate\|cut` | 構造変換 | — | ❌ | 未対応 |
| `mater` / `material` | 材料パラメータ設定 | `region material=` | 🟡 | 既知材料(si/oxide/nitride/poly/gas)のみ |

## 2. プロセスステップ

| Sprocess | 機能 | Cprocess 対応 | 状態 | 備考 |
|---|---|---|---|---|
| `implant <species> dose= energy= tilt= rotation=` | イオン注入(解析) | `implant species= dose= energy= tilt= rotation=` | ✅ | 8種(B/P/As/Sb/In/C/F/Ge)対応。元素名→記号変換 |
| `implant BF2 dose= energy=` | BF2+ 分子注入 | `implant species=B energy=E*0.226` | 🟡 | 質量比 (M_B/M_BF2≈0.226) によるエネルギースケール近似で B に分解。ドーズは B 原子数としてそのまま引継ぎ。解離・電荷分配は無視（警告コメント出力） |
| `implant ... sentaurus.mc \| crystaltrim` | MC 注入 | `implant ... method=mc` | ✅ | チャネリング・ダメージ対応 |
| `diffuse temperature= time=` (雰囲気なし) | アニール/拡散 | `diffuse temp= time=` | ✅ | ChargedFermi 相当 |
| `diffuse ... O2` | 乾式熱酸化 | `oxidize time= temp= ambient=dry` | ✅ | Deal-Grove。time/temp を引継ぎ。以前の「不活性アニールへ劣化」は解消 |
| `diffuse ... H2O\|steam\|wet` | 湿式熱酸化 | `oxidize time= temp= ambient=wet` | ✅ | 同上、wet 判定 |
| `diffuse ... O2\|H2O lateral=1` (簡易記法) | LOCOS 酸化 | `oxidize2d time= temp= ambient=` | 🟡 | 変換器独自の `lateral=` フラグで `oxidize`/`oxidize2d` を切替（Sentaurus 標準構文ではない簡易対応。実際の LOCOS 判定は nitride マスクの有無に依存するため要手動確認） |
| `deposit material= thickness=\|rate= time= type=` | デポジション | `deposit material= thickness=` | ✅ | `thickness=` はそのまま、`rate×time` は um/min×min で厚さ換算。`type=`(異方性/等方性)は現行デッキに引数が無いため反映されない(注記) |
| `etch material= thickness=\|rate= time= type=` | エッチング | `etch depth= material=` | ✅ | 同上の厚さ換算。深さのみの幾何エッチ（`etch_rate` レベルセットは変換器未対応、Python 直書き推奨） |
| `strip <material>` | 層除去 | `strip` | ✅ | レジストのみ除去（Cprocess の `strip` 実装に準拠） |
| `mask name= left= right=` | マスク定義 | implant の `x1=x2=` 窓に転用 | 🟡 | 後続 implant の窓へ反映可 |
| `photo mask= thickness=` | レジスト塗布 | `photo resist=` | ✅ | thickness をそのまま引継ぎ |
| `temp_ramp name= breakpoints=t:T,...` | 温度傾斜 | `diffuse ramp=t0:T0,t1:T1,...` | ✅ | 変換器独自の `breakpoints=`/`points=` 簡易記法のみ対応。`name=` で定義して別行の `diffuse temp_ramp=<name>` から参照する完全な Sentaurus 構文は未対応（警告） |
| `anneal` | アニール(diffuse 別名) | `diffuse` | ✅ | diffuse と同義に翻訳 |
| `silicide metal= time= temperature=` | シリサイド化 | `silicide metal= time= temp=` | ✅ | Ni/Ti 対応 |
| `stress` | 応力解析(独立コマンド) | `mechanics temp= time=` | 🟡 | `stress` という独立コマンド名は変換器未対応。`mechanics` デッキコマンドを直接使うこと |
| `epitaxy` | エピタキシャル成長 | `epitaxy thickness= temp= time=` | ❌ | 変換器のコマンド解析は未対応（Sentaurus 側構文が `epi`/`deposit type=epitaxial` など複数流儀のため）。デッキ/Python では利用可能 |
| `contact name= ...` | 電極定義 | — | ⚠️ | デバイス解析用、プロセスでは無視 |

## 3. パラメータ・制御

| Sprocess | 機能 | Cprocess 対応 | 状態 | 備考 |
|---|---|---|---|---|
| `pdbSet <Reg> <Sp> <Param> <val>` | パラメータDB上書き | `pdbset key= value=` | 🟡 | 既知キー対訳表（下記）のみ翻訳。表に無いキーは `# [warn] ... dropped` を出力（黙って消さない） |
| `pdbGet` / `pdbDelayDouble` | パラメータ取得 | — | ⚠️ | 変換器は無視（プロセス実行に影響しないため） |
| `math numThreads= numThreadsMC=` | 並列数 | implant の `threads=` | 🟡 | MC スレッド数のみ転用（未実装、無視のまま） |
| `set <var> <val>` (Tcl) | 変数定義 | 変換器内で展開 | 🟡 | `$var` を字句置換 |
| `if/else/foreach` (Tcl) | 制御構文 | — | ❌ | 変換器は非対応（警告） |
| `fset` / `define` | マクロ | 字句置換 | 🟡 | 単純定義のみ |
| `exit` / `quit` | 終了 | `stop` | ✅ | |

### pdbSet キー対訳表 (`_PDB_KEY_MAP`)

| Sentaurus 慣用名 (大小文字/`.`/`_`無視) | Cprocess ParamDB キー |
|---|---|
| `ted.frenkelSurvival` | `ted.frenkel_survival` |
| `ted.kci` | `ted.k_ci` |
| `oed.theta` | `oed.theta` |
| `oed.psiCap` | `oed.psi_cap` |
| `ox2d.cgas` | `ox2d.cgas` |
| `ox2d.nitrideLeak` | `ox2d.nitride_leak` |
| `ox2d.seedOx` | `ox2d.seed_ox_um` |
| `sige.couple` | `sige.couple` |
| `sige.dEgCoef` | `sige.dEg_coef` |
| `sige.eps0Coef` | `sige.eps0_coef` |
| `sper.actFactor` | `sper.act_factor` |
| `sper.amorphDensity` | `sper.amorph_density` |
| `sper.ea` / `sper.v0` | `sper.ea` / `sper.v0` |
| `stress.couple` | `stress.couple` |
| `stress.vr` | `stress.vr` |
| `stress.vact.<mat>` | `stress.vact.<mat>` |
| `mech.E.<mat>` / `youngsModulus.<mat>` | `mech.E.<mat>` |
| `mech.nu.<mat>` / `poissonRatio.<mat>` | `mech.nu.<mat>` |
| `mech.alpha.<mat>` | `mech.alpha.<mat>` |
| `mech.sigma0.<mat>` / `intrinsicStress.<mat>` | `mech.sigma0.<mat>` |
| `mech.tau.<mat>` | `mech.tau.<mat>` |
| `<Sym>.fi`（例 `boron.fi`） | `<Sym>.fi` |
| それ以外 | 未対応（警告コメントのみ、翻訳されない） |

Sentaurus の正式 PDB キー名は非公開のため、これは慣用名からの近似マッピングであり
網羅的ではない。実データでのキー名確認が今後必要。

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
| (Python専用) `save_state`/`load_state`/`export_device` | 状態保存/デバイス出力 | `save_state`/`load_state`/`export_device` | ✅ | Sentaurus に直接の対応コマンド無し。デッキでは利用可能だが変換器コマンドとしては未マップ（対応する `.cmd` 構文が無いため） |

## 5. 種名・単位の対応

**元素名 → Cprocess 記号**（変換器が自動、8種）:
`Boron→B, Phosphorus→P, Arsenic→As, Antimony→Sb, Indium→In, Carbon→C,
Fluorine→F, Germanium→Ge`。
**BF2 → B**（質量比エネルギースケール近似、`sprocess.py` の `_BF2_ENERGY_SCALE`
参照。ドーズはそのまま B の原子数として引継ぎ）。
（上記以外は警告）。

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

`python/cprocess/sprocess.py` を参照。ゴールデンテストは
`python/test_sprocess_converter.py`。

```python
import cprocess as cp
from cprocess.sprocess import translate_file, run_sprocess

# (1) Sprocess .cmd を Cprocess デッキ文字列へ変換
deck = translate_file("nmos.cmd")
print(deck)            # 翻訳結果（未対応行は # [unsupported] コメント、
                        # 対応外パラメータは # [warn] コメント）

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
- `region`/`line` から推定される 3D box メッシュは Sentaurus の実際の
  構造（TDR ベースの任意形状）とは一致しない（近似メッシュでの単体検証用）。
- Tcl 制御構文（`if/foreach/expr`）は非対応。事前に展開された .cmd を想定。
- `temp_ramp`/LOCOS酸化(`lateral=`)は変換器独自の簡易記法のみ対応（実際の
  Sentaurus 構文とは異なる）。
- `pdbSet` は対訳表にあるキーのみ翻訳。表に無いキーは翻訳されない
  （黙って消さず `# [warn]` を出力）。
- BF2 のエネルギースケーリングは質量比による一次近似であり、実測プロファイル
  との厳密な一致は保証しない。

---

## 参照
- Sentaurus Process チュートリアル 1〜3, 15
  <https://ghzphy.github.io/Sentaurus_Training/sp/sp_1.html>,
  <https://ghzphy.github.io/Sentaurus_Training/sp/sp_2.html>,
  <https://ghzphy.github.io/Sentaurus_Training/sp/sp_3.html>
- Sentaurus Process User Guide (N-2017.09 公開版)
- 関連: `docs/sentaurus_reference.md`, `docs/sentaurus_gap_analysis.md`,
  `docs/tasks/W1_sprocess_converter.md`
