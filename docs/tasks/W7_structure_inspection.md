# W-7: 途中構造の検査機能(レジストスタックの可視化・保存緩和)

## 目的

photo/mask_polygon の結果(レジスト形状)は隠しサイドメッシュ
`SimState::stack` + `stack_cell_mat` にのみ存在し、`save()` は本体メッシュ
のみ出力、`save_state()`/`export_device()` はレジスト存在時に throw する。
**エンジニアがマスク形状を確認する手段が皆無**(根本分析:
`docs/structure_model_root_cause.md`)。本タスクで (a) スタックの VTU 出力
`proc::save_stack`、(b) `save()` のサイドカー併記、(c) throw の警告緩和、
(d) Python 検査 API `resist_mask()` を実装する。パリティチェック
[W-7/A-7]「resist visible in saved VTU」/[W-7]「save_state with resist
stack」の 2 件を PASS に反転し、通常スイートへ移設する。

### エンジニアはこの工程の結果をどう確認するか

- **ParaView**: `sim.save("out.vtu")` がレジスト存在時に自動で
  `out_stack.vtu` を併記。セルデータ `Material_si0_resist1_open2`
  (0=Si 基板, 1=レジスト, 2=現像開口)で形状を色分け表示できる。
  配列名自体が材料マッピングの説明になっている。
- **matplotlib / numpy**: `sim.resist_mask()` が (N,4) numpy 配列
  (セル中心 x,y,z [µm] + 材料インデックス)を返す。
  `plt.scatter(a[:,0], a[:,1], c=a[:,3])` で即座に平面マスクを確認。
- **ログ**: save/save_stack/save_state がレジストセル数・開口セル数・
  出力パス・警告をログ文字列に出す。

## 現状コード

- `include/cprocess/deck.hpp` — `SimState`: `has_stack`, `stack` (Mesh),
  `stack_cell_mat` (0=Si/1=resist/2=open), `stack_resist_z0`
- `src/process.cpp` — `photo()`/`mask_polygon()` がスタックを構築・現像。
  `save()` は本体メッシュのみ `write_vtu`(int 配列 "Region" の書式が流用元)。
  `save_state()`/`export_device()` は `has_stack` で throw
- `src/vtk_writer.cpp` — `write_vtu(path, mesh, scalars, int_scalars,
  point_scalars={})`(P3-h で point_scalars 追加。今回は int_scalars のみ使用)
- `tests/test_sprocess_parity.cpp` — 対象 2 チェック(WILL_FAIL 運用)

## 実装手順

1. `proc::save_stack(st, path, log)`:
   `!has_stack` なら throw。`write_vtu(path, st.stack, {},
   {{"Material_si0_resist1_open2", &mat}})` でスタックメッシュを出力。
   配列名に "resist" を含める(マッピングの自己記述 + 検出基準を正直に充足:
   ファイルにはレジストセル形状そのものが入る)。
2. `proc::save(st, path, log, include_stack=true)`:
   **既定でサイドカー併記**。`include_stack && has_stack` のとき
   `<path から .vtu を除いた名前>_stack.vtu` を save_stack で併記しログ。
   - **設計判断**: 計画草案は「既定 off のオプション」だったが、パリティ
     チェックは既定の `save()` 呼び出しでレジスト形状が見えることを要求する。
     レジスト存在下で構造を保存するエンジニアはレジストを見たいはず、が
     本タスクの発端(検査手段の欠如)なので、**既定 on** を採用。
     `include_stack=false` で従来出力(バイト同一)に戻せる。
3. `save_state()`/`export_device()` の throw 緩和:
   **throw をやめ、警告ログ + レジスト抜きで保存続行**に変更。
   - **設計判断**: 計画 W-7 (4) は「既定は throw、`force=` で緩和」だったが、
     パリティチェックは素の呼び出しが throw しないことをアサートしており
     矛盾する。photo/mask 状態は安価に再構築できるためデータ損失は軽微で、
     silent-throw こそが使い勝手の苦情の原因だったので、**警告 + 続行**を
     採用し計画本文を修正した(`force=` パラメータは追加しない)。
     警告文にスタックが永続化されない旨を明記する。
4. Python 三点セット(CLAUDE.md 必須):
   - pybind `proc_save_stack` / `proc_resist_mask`、`proc_save` に
     `include_stack` 引数追加
   - `Simulation.save(path, include_stack=True)` / `save_stack(path)` /
     `resist_mask()`((N,4) numpy 配列、中心座標 µm + 材料インデックス。
     クエリのため self でなくデータを返す)
5. パリティチェック 2 件を移設(WILL_FAIL 運用):
   - resist 可視化 → `tests/test_photo.cpp`(save の既定サイドカー +
     save_stack の内容を検証)
   - save_state 非 throw → `tests/test_state_io.cpp`(警告ログ +
     ラウンドトリップ検証)
   - `docs/tasks/README.md` のパリティ表を 修正済/移設 に更新、
     サマリ表記 10 → 8。`docs/IMPLEMENTATION_PLAN_v2.md` W-7 (4) を修正。

## テスト仕様

- C++ `tests/test_photo.cpp`: photo(0.2µm)+mask_polygon(半面開口)後、
  `save()` 既定呼び出しで `<base>_stack.vtu` が生成され、"resist" 文字列を
  含み、セル数がベースメッシュ超であること。`save_stack` 単体も同内容。
  スタックなしでの `save_stack` は throw。`include_stack=false` では
  サイドカーを生成しない。
- C++ `tests/test_state_io.cpp`: レジスト存在下の `save_state` が throw
  せず警告をログし、`load_state` 後のベース状態が一致すること。
- Python `python/test_comprehensive.py`: `save_stack`(ファイル生成 +
  "resist" 含有)、`resist_mask()`(形状 (N,4)、材料 {0,1,2}、開口セルの
  x 座標がマスク窓内)、`save(include_stack=...)` の挙動。

## DoD

- パリティ 2 チェックが移設先で PASS、`sprocess_parity` は 0/8 で
  WILL_FAIL 維持(ctest green)。
- C++ 全 52 テスト PASS、`python3 python/test_comprehensive.py` /
  `test_simulation.py` PASS。
- 複雑ポリゴンマスクを photo → mask_polygon → save_stack で ParaView
  確認できる(golden flow / test_photo で担保)。

## やらないこと

- スタックの CPRC1 シリアライズ(save_state はレジスト抜き保存 + 警告。
  必要になれば別タスク)。
- 構造モデル不変条件(CLAUDE.md)の根本解消 — スタックを本体表現へ統合
  する作業は A-7 の範囲。本タスクは検査手段の提供のみ。
- GDS 読込マスクの追加テスト(既存 test_gds でカバー)。
