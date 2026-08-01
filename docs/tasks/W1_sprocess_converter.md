# W-1: sprocess.py 変換器の全面更新

## 目的

`python/cprocess/sprocess.py`（Sentaurus Process `.cmd` → Cprocess デッキ変換器）
は 2026-06-15 時点のエンジン能力（`mesh/region/init/implant/diffuse/mask/struct/set`
のみ）を前提に書かれている。以降のスプリントで `oxidize`/`oxidize2d`/`deposit`/
`etch`/`silicide`/`epitaxy`/`sper`/`mechanics`/`refine`/`pdbset`/`save_state`/
`load_state`/`export_device` がデッキコマンドとして実装され（`src/deck.cpp`
2026-06-15 以降）、Python `Simulation` も 37+ メソッドまで拡張されたが、
変換器は追随しておらず、実装済み機能を次のように劣化・無視していた:

- `diffuse ... O2|H2O` → 酸化を無視した不活性アニールに劣化
- `deposit`/`etch`/`silicide`/`temp_ramp` → unsupported（未実装扱い）
- `photo`/`strip` → ignore（構造変化を捨てる）
- `pdbSet` → 警告のみで完全ドロップ
- 対応ドーパント種が B/P/As/Sb の 4 種のみ（`src/materials.cpp` の
  `kDopants` は B/P/As/Sb/In/C/F/Ge の 8 種）

このタスクは変換器を現行エンジン実装に追随させ、「エンジンに実装済みの
コマンドを 1 つも unsupported と誤判定しない」状態にする。

## 現状コード

- `python/cprocess/sprocess.py` — 変換器本体。`_cmd_*` 関数群 + `_IGNORE`/
  `_HARD` セット方式のディスパッチ。
- `src/deck.cpp` — 現行デッキコマンド（`run_deck` の if/else 連鎖）:
  `mesh, region, init, implant, photo, mask, strip, bc, diffuse, oxidize,
  epitaxy, deposit, silicide, etch, pdbset, oxidize2d, sper, mechanics,
  refine, save_state, load_state, export_device, save, print, stop`
  （24 コマンド。`docs/IMPLEMENTATION_PLAN_v2.md` の「16 コマンド」という
  記述は W-3 で増設済みのため古い）。
- `python/cprocess/simulation.py` — Python `Simulation` API（変換後の
  デッキ文字列が最終的に呼び出す先の等価物）。
- `src/materials.cpp` の `kDopants`: `boron/B, phosphorus/P, arsenic/As,
  antimony/Sb, indium/In, carbon/C, fluorine/F, germanium/Ge` の8種。
- ParamDB キー一覧（`grep -rn 'ParamDB::instance().get(\|db.get(' src/*.cpp`
  で収集）: `oed.theta, oed.psi_cap, ox2d.cgas, ox2d.nitride_leak,
  ox2d.seed_ox_um, sige.couple, sige.dEg_coef, sige.eps0_coef,
  sper.act_factor, sper.amorph_density, sper.ea, sper.v0, stress.couple,
  stress.vr, stress.vact.<mat>, ted.frenkel_survival, ted.k_ci,
  pd.*（点欠陥）, cl.b.*/cl.as.*（クラスタ化）, mech.E./mech.nu./
  mech.alpha./mech.sigma0./mech.tau.<mat>, I.*, <Sym>.fi` など。
  Sentaurus 側のキー名は非公開のため、既知の慣用名（TCAD 文献・
  Sentaurus User Guide 目次に現れる代表的な PDB パス）から Cprocess の
  対応キーへの**近似マッピング**を用意し、対応外は警告のまま残す。

## 実装方針

1. **`diffuse ... O2|H2O`** — 従来の「警告して不活性アニール」をやめ、
   雰囲気を検出したら `oxidize`（`oxidize2d` は明示 `lateral=`/2D 指定が
   あるときのみ；既定は 1D `oxidize`）へ変換する。`O2` → dry、`H2O`/
   `steam`/`wet` → wet と判定し、`time=`/`temperature=` をそのまま
   引き継ぐ。雰囲気なしの `diffuse` は従来通り `diffuse` のまま。
2. **`deposit`/`etch`** — `thickness=` があれば直接 `deposit
   material= thickness=`／`etch depth= material=` に変換。`rate=` +
   `time=` の組み合わせは `thickness = rate × time` を変換器内で計算
   （µm/min × min = µm）。`type=`（anisotropic/isotropic）は現行デッキに
   引数が無いため注記コメントを残し警告はしない（構造上の近似）。
3. **`silicide`** — `metal=`/`temperature=`/`time=` を `silicide metal=
   temp= time=` に変換。
4. **`photo`/`strip`** — `photo mask= thickness=` を `photo resist=` へ、
   `strip` をそのまま `strip` へ変換（従来の ignore をやめる）。
5. **`temp_ramp`** — 複数ブレークポイントを `diffuse ramp=t0:T0,t1:T1,...`
   （`src/deck.cpp` の `parse_ramp` 形式、time は累積分、`ramp[0]` は
   t=0）に変換。
6. **`pdbSet`** — `pdbSet <Reg> <Sp> <Param> <val>` の `<Param>` 部分を
   既知対訳表（下記）で Cprocess `ParamDB` キーへマップし `pdbset key=
   value=` を出力。対訳表に無いキーは翻訳を諦め、明示的な警告コメント
   （`# [warn] pdbSet '<name>' has no known Cprocess ParamDB equivalent —
   dropped`）を出す（黙って消さない）。
7. **8 種ドーパント対応** — `_SPECIES` を B/P/As/Sb/In/C/F/Ge の8種に拡張。
   `BF2`/`bf2` は Cprocess に BF2 という核種が無いため、**質量ベースの
   分解近似**として B にマップし、実効打ち込みエネルギーを
   `E_B = E_BF2 * (M_B / M_BF2)`（一次近似、SRIM 等で使われる質量比
   スケーリングと同じ考え方；電荷状態・分子解離の効果は無視する近似で
   あることをコメントで明記）で換算する。`M_BF2 ≈ 48.6 amu`、
   `M_B ≈ 11.0 amu` を使用（比 ≈ 0.226）。ドーズはそのまま B の原子数
   として引き継ぐ（BF2+ 1 個 = B 原子 1 個という実務上の慣用に合わせる）。

## パラメータDB対訳表（pdbSet）

変換器内 `_PDB_KEY_MAP` に実装。Sentaurus 側の慣用パラメータ名 → Cprocess
ParamDB キー:

| Sentaurus 慣用名 (小文字比較) | Cprocess キー | 備考 |
|---|---|---|
| `oxide.dry.rate` / `oxidedryrate` 系 | 未対応 | Deal-Grove 係数はハードコード。警告 |
| `frenkelsurvival` / `ted.frenkelsurvival` | `ted.frenkel_survival` | |
| `ted.kci` / `interstitial.kci` | `ted.k_ci` | |
| `oed.theta` / `oxidant.theta` | `oed.theta` | |
| `oed.psicap` / `oed.psi_cap` | `oed.psi_cap` | |
| `ox2d.cgas` / `oxidant.cgas` | `ox2d.cgas` | |
| `ox2d.nitrideleak` | `ox2d.nitride_leak` | |
| `ox2d.seedox` | `ox2d.seed_ox_um` | |
| `sige.couple` | `sige.couple` | |
| `sige.dEgcoef` / `sige.bandgapnarrowing` | `sige.dEg_coef` | |
| `sige.eps0coef` | `sige.eps0_coef` | |
| `sper.actfactor` | `sper.act_factor` | |
| `sper.amorphdensity` | `sper.amorph_density` | |
| `sper.ea` | `sper.ea` | |
| `sper.v0` | `sper.v0` | |
| `stress.couple` | `stress.couple` | |
| `stress.vr` | `stress.vr` | |
| `stress.vact.<mat>` | `stress.vact.<mat>` | material 名は透過 |
| `mech.e.<mat>` / `youngsmodulus.<mat>` | `mech.E.<mat>` | |
| `mech.nu.<mat>` / `poissonratio.<mat>` | `mech.nu.<mat>` | |
| `mech.alpha.<mat>` | `mech.alpha.<mat>` | |
| `mech.sigma0.<mat>` / `intrinsicstress.<mat>` | `mech.sigma0.<mat>` | |
| `mech.tau.<mat>` | `mech.tau.<mat>` | |
| `<Sym>.fi` (例 `boron.fi`) | `<Sym>.fi` | 種記号は `_SPECIES` で正規化 |
| それ以外 | — | `# [warn] ... dropped` |

## テスト仕様（`python/test_sprocess_converter.py`）

ゴールデンテスト5本 + 単体アサーション:

1. **LOCOS + well + S/D フロー** — `line`/`region`/`init`/`mask`/
   `implant`（well）/`diffuse O2`（フィールド酸化）/`implant`（S/D，BF2）/
   `diffuse` を含む代表レシピ → 変換後デッキに `oxidize`、
   `implant species=B`（BF2分解）が含まれ、`# [unsupported]` 行が
   0 件であることを確認。
2. **`diffuse` + `O2`（dry酸化）** → 出力に `oxidize time=... temp=...`
   が含まれ `ambient` を暗黙 dry（`wet` 指定なし）扱いで正しく生成する
   こと。旧来の「不活性アニール」警告コメントが**出ない**こと。
3. **`diffuse` + `H2O`（wet酸化）** → 出力に `oxidize ... ambient=wet`相当
   （`cp.Simulation.oxidize(wet=True)` に対応する `wet=` 引数、または
   `src/deck.cpp` の `ambient=wet` キーワード）が含まれること。
4. **`deposit` + `etch` ペア**（`rate=`×`time=` 指定含む）→
   `deposit material=... thickness=...` と `etch depth=... material=...`
   が両方出力され、`thickness = rate*time` の数値が正しく換算される
   こと。
5. **`pdbSet` オーバーライド** — 既知キー（例 `sper.actfactor`）は
   `pdbset key=sper.act_factor value=...` に翻訳され、未知キーは
   `# [warn] ... dropped` コメントになること（両方を1レシピでテスト）。

各テストは `translate()` の出力が「構文的に有効なデッキ」であることも
検証する（各 `deposit`/`etch`/`implant`/`diffuse`/`oxidize`/`pdbset` 行が
`src/deck.cpp` の `Cmd` パーサ規約 `key=value` トークン形式に従っている
ことを正規表現で確認）。可能なら `cp.run_deck()` に実際に通して例外なく
実行できることも確認する（C++ 拡張がビルド済みの環境でのみ、`try/except
ImportError` でスキップ可）。

## 完了条件 (DoD)

- `docs/sprocess_command_map.md` が新しい対訳表と一致している。
- 変換器がエンジン実装済みコマンド（`oxidize/oxidize2d/deposit/etch/
  silicide/epitaxy/sper/mechanics/refine/pdbset/save_state/load_state/
  export_device`）を1つも `# [unsupported]` にしない。
- 代表 LOCOS+well+S/D レシピが警告最小（意図しない unsupported 0件）で
  変換される。
- 8 種ドーパント + BF2 分解に対応。
- `python/test_sprocess_converter.py` の5ゴールデンテストが pass。
- `python3 python/test_comprehensive.py && python3 python/test_simulation.py`
  にリグレッションが無い。

## エンジニアによる結果の検査方法（必須項目）

変換結果を検査する主な方法は3通り:

1. **デッキ文字列の目視/diff**: `translate_file("foo.cmd")` の戻り値を
   そのまま `print()` するか `-o out.deck` でファイル出力し、
   `# [unsupported]` / `# [warn]` 行の有無を `grep` で確認する（CI でも
   `assert "[unsupported]" not in deck` の形で機械チェック可能）。
2. **実行して構造を確認**: `run_sprocess("foo.cmd")` で実際にエンジンを
   走らせ、返る `SimState` を `sim.save("out.vtu")` 相当（`cp.run_deck`
   の `save` コマンド経由）で ParaView 出力し、酸化膜厚・エッチ深さ・
   ドーパントプロファイルが Sentaurus 側の期待値と定性的に一致するかを
   目視する。
3. **ユニットテスト**: `python/test_sprocess_converter.py` のゴールデン
   テストが「入力 `.cmd` 断片 → 期待デッキ行を含む」ことを CI で継続的に
   保証する。将来 `src/deck.cpp` にコマンドが追加/変更された場合、この
   テストが真っ先に失敗することで変換器の追随漏れを検出できる。

## やらないこと

- Tcl 制御構文（`if/foreach/expr`）の実行時展開（従来通り非対応、
  事前展開された `.cmd` を前提とする）。
- `pdbSet` の全キー網羅（対応表にないキーは警告のまま。Sentaurus の
  正式 PDB キー名は非公開のため実データでの追加検証が今後必要）。
- `docs/sentaurus_gap_analysis.md` 等、`sprocess_command_map.md` 以外の
  広範な文書同期（W-2 の範囲）。
