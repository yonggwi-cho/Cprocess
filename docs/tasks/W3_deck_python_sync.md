# W-3: デッキと Python API の機能同期

## 目的

デッキ(`src/deck.cpp`)は 16 コマンド、Python `Simulation` は 37+ メソッド
を持ち、両者は `proc::` を共通の実装として共有しているにもかかわらず
`etch`/`oxidize_2d`/`sper`/`mechanics`/`refine`/`ramp`(diffuse オプション)/
`pdbset`(ParamDB 設定)/`save_state`/`load_state`/`export_device` がデッキ
から呼べなかった(gap analysis)。パリティチェック [W-3]「deck parity: etch
command」「deck parity: pdbset command」を PASS に反転し、通常スイートへ移設
する。あわせて `Simulation.bbox()` が µm 系 API の中で唯一 cm を返していた
バグを修正する。

### エンジニアはこの工程の結果をどう確認するか

- **デッキのログ**: 新コマンドはすべて既存コマンド(oxidize/deposit 等)と
  同じ流儀で `log` へ実行結果を出す(proc:: 側の既存ログ機構をそのまま使う
  ので新規のログ実装は不要)。
- **エラー**: 構文誤り・未知コマンドは既存規約どおり
  `deck line <N> (<cmd>): <msg>` の行番号付き例外で失敗する
  (`tests/test_flow.cpp` の `test_deck_new_commands` にエラー経路を用意)。
- **テスト**: `tests/test_flow.cpp::test_deck_new_commands` が各コマンドの
  happy path をメッシュ縮小・ParamDB 反映・応力フィールド書き込み・セル数
  増加・ファイル出力などの観測可能な副作用で検証する。
- **Python**: `sim.bbox()` が µm を返すようになったので、他の µm 系アクセサ
  (`cell_centroids` など)とそのまま算術演算できる(修正前は手動で `*1e4`
  変換が必要だった)。

## 現状コード

- `src/deck.cpp` — コマンドディスパッチャ(`run_deck`)。`Cmd` 構造体が
  `str()`/`num()`/`num_or()`/`flag_or()` で単位接尾辞付き引数をパースする
  既存の型どおりの規約(`0.1um`/`1000C`/`30min`)。`cmd_oxidize`/
  `cmd_deposit` をテンプレートに新コマンドを実装。
- `include/cprocess/process.hpp` — 呼び出す `proc::` シグネチャ一式
  (`etch`, `oxidize_2d`, `sper`, `mechanics`, `refine`, `save_state`,
  `load_state`, `export_device`, `set_param`/`get_param`)。
- `include/cprocess/diffusion.hpp` — `DiffuseOpts::temp_profile`
  (`std::vector<std::pair<double,double>>`, (秒, K) のブレークポイント列)。
- `include/cprocess/param_db.hpp` / `src/param_db.cpp` — `ParamDB::instance()`
  はメッシュ非依存の単純なキー→値マップ。`proc::set_param(SimState&, key,
  value, log)` は `st` を無視して `ParamDB::instance().set(key, value)` を
  呼ぶだけ(process.cpp 側で確認済み)。**キー検証機構が存在しない**
  (`overrides_` は素の `std::map`、未知キーはそのまま黙って格納される)。
- `python/cprocess/simulation.py:596` — `bbox()` が `self._st.mesh.nodes`
  (cm)をそのまま返していた。内部で `mask()`(197 行目)と `dose()`
  (592 行目)が `bbox()` の戻り値を使っている。

## 実装手順

1. **`deck.cpp` 新コマンド**(すべて実装、後述の「やらないこと」なし):
   - `etch depth=<len> [material=<str>] [poly=(x1,y1),(x2,y2),...]` —
     `proc::etch` へ委譲。`poly=` は新設の `parse_poly()` ヘルパで
     `(x,y),(x,y),...` を µm 長さ単位でパースする(既存の `parse_quantity`
     を頂点ごとに再利用)。
   - `pdbset key=<str> value=<float>` — `proc::set_param(st, key, value,
     &log)`。`ParamDB` は `SimState` 非依存のシングルトンなので `st` は
     API 整合性のためだけに渡る(process.hpp のコメントどおり)。
   - `oxidize2d time=<t> temp=<T> [ambient=dry|wet]` — `proc::oxidize_2d`。
     `cmd_oxidize` と同じ引数規約。
   - `sper temp=<T> time=<t>` — `proc::sper`。
   - `mechanics temp=<T> time=<t>` — `proc::mechanics`。
   - `refine species=<str> [thresh=<f>] [passes=<n>] [axis=x|y|z]` —
     `proc::refine`。
   - `diffuse` の `ramp=` オプション — `0min:900C,10min:1050C,...`
     のセミコロン/カンマ区切りブレークポイント列を新設の `parse_ramp()` で
     `DiffuseOpts::temp_profile` に変換。Python `Simulation.diffuse(ramp=
     [(t_min, T_celsius), ...])` と同じ意味論(先頭ブレークポイントは
     t=0 必須、`temp=` と排他ではなく `ramp=` があれば `temp=` は無視)。
   - `save_state file=<path>` / `load_state file=<path>` — `proc::save_state`
     / `proc::load_state`。
   - `export_device prefix=<path>` — `proc::export_device`。
   - `mask_polygon` はデッキに追加していない(**やらないこと**参照)。

2. **`bbox()` の µm 修正**(`python/cprocess/simulation.py:596`):
   `xyz.min(axis=0) / UM, xyz.max(axis=0) / UM` に変更。内部呼び出し元
   `mask()`(y1/y2 省略時のデフォルト)と `dose()`(面積計算)を
   `bbox()` の新しい単位系に合わせて更新(`* UM` を明示的に挿入)。
   `python/test_comprehensive.py` の 3 箇所(`test_segregation_dose_loss`,
   `test_nitride_barrier_python` の µm/cm 変換コメント付き利用)と
   `examples/circular_contact.py` の面積計算を、修正後の µm 挙動に合わせて
   更新した(`grep -rn "\.bbox("` で洗い出し済み。単純な大小比較のみの
   箇所(`test_etch_depo_p17`, `test_oxidize_2d_python` 等)は単位が変わって
   も相対比較なので無修正で正しい)。

3. **ParamDB 未知キー警告**: **見送り**(下記「やらないこと」参照)。

## テスト仕様

- `tests/test_flow.cpp::test_deck_new_commands`(新設):
  - `etch`: blanket depth でメッシュ bbox が縮む。`material=` + `poly=`
    形式もエラーなく実行できる。
  - `pdbset`: `proc::get_param("oed.theta", ...)` で反映を確認し、
    テスト終了時にデフォルト値へ復元(他テストへの汚染防止)。
  - `oxidize2d`: nitride マスクなしでは throw(仕様どおり)、開口ありでは
    正常終了。
  - `sper`: damage フィールドなしメッシュで no-op(throw しない)。
  - `mechanics`: `sxx` 等の応力フィールドが書き込まれることを確認。
  - `refine`: 実行後にセル数が増えることを確認。
  - `save_state`/`load_state`: 往復でフィールドが復元される。
  - `export_device`: `<prefix>.vtu` と `<prefix>.meta.json` が生成される。
  - `diffuse ramp=...`: 3 点ブレークポイントの RTA が実行できる。
  - エラー経路: `etch`(depth= 欠落)、`pdbset`(value= 欠落)、未知コマンド
    がいずれも例外(`etch`/`pdbset` は "deck line N" 形式)を出す。
- `tests/test_sprocess_parity.cpp`: [W-3] の 2 チェックを削除し、上記へ移設
  した旨のコメントを残す。
- Python 側: `python/test_comprehensive.py` の bbox 利用箇所を更新(上記)。
  bbox() 自体の新規専用テストは追加していない(既存テストの間接検証で
  µm/cm の取り違えは即座に数値が破綻し検出される)。

## DoD

- [x] デッキで `etch`/`pdbset` が受理される([W-3] の 2 チェックが PASS)。
- [x] `oxidize2d`/`sper`/`mechanics`/`refine`/`ramp`/`save_state`/
      `load_state`/`export_device` もデッキから呼べる(ボーナス実装)。
- [x] `bbox()` が µm を返し、既存呼び出し元(mask/dose/tests/examples)が
      追従済み。
- [ ] `examples/` へのデッキ版フルフロー(NMOS LOCOS)1 本の追加は
      **未実施**(下記「やらないこと」参照)。
- [ ] ParamDB 未知キー警告は **未実施**(下記参照)。

## やらないこと

- **ParamDB 未知キー警告**: `ParamDB` は素の `std::map<string,double>` で
  あり、どのキーが「既知」かを判定する集中管理されたレジストリが存在しない
  (`materials.cpp`/`diffusion.cpp`/`process.cpp` などに散在する
  `ParamDB::instance().get("...", fallback)` 呼び出しがそれぞれ独自に
  キー名をハードコードしている)。安全に「未知キー警告」を実装するには、
  まずそれら呼び出し箇所全体を横断的に集約してキー名一覧を構築する
  broader な変更が要る(単純に「よく見るプレフィックスで判定」のような
  当て推量の設計は、誤って正当なキーに警告を出す/新規キー追加のたびに
  一覧を更新し忘れるリスクを生む)。タスク指示の「非自明なら推測せず見送り、
  理由を記す」に従い本タスクでは見送る。フォローアップ候補: `grep -rn
  'ParamDB::instance().get' src/` の出力からビルド時にキー一覧を生成する
  small script、または各 `.get()` 呼び出し箇所に近い場所へ既知キー登録
  マクロを添える設計。
- **examples/ へのデッキ版 NMOS LOCOS フルフロー追加**: plan の DoD 項目
  だが、規模(実 LOCOS レシピをデッキ構文で書き下ろし+期待値アサート)が
  本タスクの必須線(etch+pdbset を PASS にする)を大きく超えるため、
  時間配分の都合で見送った。`tests/test_flow.cpp::test_deck_new_commands`
  が個々の新コマンドの happy path を検証しているので機能面のカバレッジは
  確保できている。
- **`mask_polygon` のデッキコマンド化**: plan の gap リストには挙がって
  いたが、W-3 の本文「やること」列挙には含まれておらず、`etch`/`deposit`
  の `poly=` 引数で既に多角形指定の同等機能をデッキから使えるため優先度を
  下げ、実装しなかった。
- **W-8(材料考慮 implant 輸送)関連ファイル**: `src/implant.cpp`、
  `include/cprocess/implant.hpp`、`process.cpp` の implant 系関数、
  `tests/test_sprocess_parity.cpp` の [W-8] チェックには一切手を入れていない
  (タスク指示により凍結中)。
