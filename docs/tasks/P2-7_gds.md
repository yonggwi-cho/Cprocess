# P2-7: GDSII レイアウト入力

## 目的

マスクポリゴンの手入力 (`mask_polygon` に座標リストを渡す) を、
実レイアウト (GDSII) からの読込に置き換えられるようにする。
外部ライブラリ依存なしの最小 GDSII バイナリリーダを実装し、
指定レイヤの BOUNDARY ポリゴンを `mask_polygon`/`etch(poly=)` に
そのまま渡せる形で返す。

## 現状コード

- `src/process.cpp` — `mask_polygon(st, poly)` / `etch(..., poly)` は
  `std::vector<std::pair<double,double>>` (cm) を受ける — 出力先はこれ
- 依存追加は不可 (CMakeLists は外部ライブラリなしが方針)

## GDSII フォーマット (実装者向け要点)

レコード列: [2B 長さ (ビッグエンディアン, レコード全体)] [1B タイプ]
[1B データ型] [ペイロード]。必要なレコードタイプ (hex):
- `0x0102` HEADER, `0x0202` BGNLIB, `0x0206` LIBNAME
- `0x0305` UNITS — ペイロードは 8B 実数×2 (GDS 独自の excess-64
  浮動小数; **変換式を実装**: 符号 1b、指数 7b (基数 16、バイアス 64)、
  仮数 56b; `value = (−1)^s · mantissa/2^56 · 16^(exp−64)`)。
  1 番目 = DB 単位のユーザ単位比 (通常 1e-3)、2 番目 = DB 単位の
  メートル値 (通常 1e-9)
- `0x0502` BGNSTR, `0x0606` STRNAME, `0x0700` ENDSTR
- `0x0800` BOUNDARY, `0x0D02` LAYER (int16), `0x0E02` DATATYPE,
  `0x1003` XY (int32 ペア列、閉多角形: 最終点=先頭点), `0x1100` ENDEL
- `0x0400` ENDLIB
その他のレコード (SREF/AREF/PATH/TEXT 等) は**スキップ** (長さ分読み飛ばす)。

## 実装手順

1. `include/cprocess/gds_reader.hpp` + `src/gds_reader.cpp` 新設:
   ```cpp
   // layer < 0 は全レイヤ。戻り値は cm 単位の閉多角形リスト
   // (末尾の重複点は取り除いて返す)。ファイル不正は runtime_error。
   std::vector<std::vector<std::pair<double,double>>>
   read_gds(const std::string& path, int layer = -1);
   ```
   - ビッグエンディアン読取ヘルパ (`u16be`/`i32be`) を実装
   - UNITS の meters-per-DB-unit を保持し、XY 座標 × (m/DBU) × 100
     で cm へ変換
   - BOUNDARY..ENDEL の間で LAYER と XY を拾い、layer フィルタ一致時に
     結果へ追加。階層 (SREF) は非対応 — 見つけたら警告カウントし
     読み飛ばす (ログは戻り値に含めないため `std::cerr` ではなく
     読み飛ばし数を返す構造体にする案もあるが、シンプルに
     ポリゴンのみ返し、SREF 数は無視する。仕様として明記)
2. proc 層 (`process.hpp`/`process.cpp`):
   ```cpp
   std::vector<std::vector<std::pair<double,double>>>
   load_gds(const std::string& path, int layer, std::ostream* log = nullptr);
   ```
   (read_gds の薄い委譲 + ログ 1 行「N polygons from layer L」)
3. pybind `proc_load_gds(state 不要 — state を取らない free 関数として
   `m.def("load_gds", ...)` で公開)。
   Python `Simulation.load_gds(path, layer=-1, scale=1.0)`:
   cm → µm 変換して `[(x,y), ...]` のリストのリストを返す
   (`scale` は追加倍率、既定 1)。**query メソッドなので self ではなく
   データを返す** (CLAUDE.md の例外規定に従う)
4. テスト用 GDS ライタ (テストローカル、`tests/test_gds.cpp` 内の
   static 関数): HEADER/BGNLIB/LIBNAME/UNITS(1e-3, 1e-9)/BGNSTR/STRNAME/
   BOUNDARY(layer2 の 1µm 角 4 点矩形)/BOUNDARY(layer5 の三角形)/
   ENDSTR/ENDLIB を正確なバイト列で書き出す。excess-64 実数の
   エンコーダも実装 (1e-3 と 1e-9 の 2 値のみ正しければよい —
   検算値をコメントに記す: 1e-9 = 0x3944B82FA09B5A54 近傍。
   実装はエンコーダを書くのではなく**この 2 定数のバイト列を
   直接埋め込む**方式でよい: 1e-3 → `0x3E 0x41 0x89 0x37 0x4B 0xC6
   0xA7 0xF0`、1e-9 → `0x39 0x44 0xB8 0x2F 0xA0 0x9B 0x5A 0x54`)
5. CMakeLists: `gds_reader.cpp` を cprocess_core に追加、
   `test_gds` を foreach に追加

## テスト仕様

`tests/test_gds.cpp`:

1. **ラウンドトリップ**: テスト内ライタで矩形 (layer 2, (0,0)-(1µm,1µm)) と
   三角形 (layer 5) を書き、`read_gds(path, 2)` → 1 ポリゴン、
   4 頂点、各座標が 1e-4 cm (=1µm) スケールで誤差 < 1e-9 cm
2. **レイヤフィルタ**: `read_gds(path, 5)` → 三角形のみ (1 ポリゴン
   3 頂点)。`read_gds(path, 1)` → 0 ポリゴン。`layer=-1` → 2 ポリゴン
3. **異常系**: 存在しないファイル → runtime_error;
   HEADER のない不正ファイル (テキストを書いたもの) → runtime_error
4. **閉多角形処理**: GDS の末尾重複点が除去されて返る (4 頂点であって
   5 頂点でない)

Python (`python/test_comprehensive.py`):
5. C++ テストが生成した .gds を tempfile に再現 (バイト列を Python 側
   にも埋め込む) → `sim.load_gds(path, 2)` → µm 座標一致;
   photo → `mask_polygon(polys[0])` → MC 注入でポリゴン外がブロック
   される (circular_contact と同じ判定、比 > 3x)

## 完了条件 (DoD)

- [ ] read_gds / proc::load_gds / pybind / Simulation.load_gds の三層完備
- [ ] 上記 C++ 4 テスト + Python テスト PASS、既存テスト回帰なし
- [ ] SREF/PATH 等はスキップされ、クラッシュしない (三角形テストの
      ファイルに PATH レコードを 1 つ混ぜて検証)
- [ ] コミットメッセージに `P2-7` を含める

## やらないこと

- SREF/AREF の階層展開・座標変換
- OASIS フォーマット
- ポリゴンのブーリアン演算 (穴あき・自己交差の正規化)
- GDS 書き出し (テストローカルのライタのみ)
