# W-2: 陳腐化文書の更新

## 1. 目的

`docs/sentaurus_gap_analysis_v2.md`・`docs/structure_model_root_cause.md`・
`docs/IMPLEMENTATION_PLAN_v2.md`・`docs/sprocess_command_map.md` は、それぞれ
作成時点(HEAD `308c2ae`/`1fe38a9` 付近)のエンジン能力を前提に書かれている。
その後 W-1/W-3/W-4/W-7/W-8/C-1/C-3 と TED 安定性修正(b7112bd)が着地し、
v2 ロードマップの第 1 群(配線・更新)と第 2 群(校正データ)の大半が解消済み
となった。文書がこれに追随していないと、読者(次の実装エージェント含む)が
「まだ無い」と誤認して重複実装したり、逆に「もう直った」ものを見落とす。
本タスクは実装済みギャップに現状注記を付け、まだ残るギャップと明確に区別する。

## 2. 現状コード(対象文書と根拠)

- `docs/sentaurus_gap_analysis_v2.md`: 1.1/1.2/1.3/2.1/2.3/3.1/3.3/§4 の各表が
  W-1〜W-4, W-7, W-8, C-1, C-3 着地前の状態を記載。
- `docs/structure_model_root_cause.md`: §2「直接原因」の記述自体は歴史的事実として
  正しいが、W-7/W-8 で是正済みである旨の注記が無い。
- `docs/IMPLEMENTATION_PLAN_v2.md`: Sprint 1/3 の各タスク節に完了マークが無い。
- `docs/sprocess_command_map.md`: W-1 (commit `79106fe`) で既に全面改訂済み
  (更新日ヘッダに明記)。本タスクでは重複改訂せず、他文書との整合のみ確認。
- `docs/tasks/README.md`: 既に各タスク着地時にステータス更新されており
  (W-1〜W-8, C-1, C-3 の状態表・sprocess_parity 表・benchmark 表が反映済み)、
  本タスク時点で追加更新の必要はほぼ無いことを確認(要 diff 確認のみ)。

## 3. 実装手順

1. `sentaurus_gap_analysis_v2.md`: 冒頭に「本書の読み方」節を追加し、
   解消済みカテゴリ(C, D の大部分 + A の一部)と残存カテゴリ(B, A の一部, C の残り)
   を一覧化。各表の該当行に `→ 解消 (W-n, commit <hash>)` 注記を追記。
   §4 ランキング表にも同様の状態列/注記を追加。
2. `structure_model_root_cause.md`: §2 冒頭に「本節は歴史的記録。W-7/W-8 で
   是正済み。現状は `docs/tasks/W7_structure_inspection.md` /
   `docs/tasks/W8_material_aware_implant.md` を参照」の注記ボックスを追加。
   §4.1 の是正タスク一覧にも完了マークを追加。
3. `IMPLEMENTATION_PLAN_v2.md`: 冒頭サマリ表の下に、W-1/W-3/W-4/W-7/W-8/C-1/C-3
   の完了状態を示す表を追加(commit ハッシュ・実装内容の一言・DoD 達成状況)。
   C-2 は着手時点の状況(TED 安定性のみ修正済み、クラスタ校正は継続中)を明記し、
   コミット直前に再確認して必要なら更新。
4. `sprocess_command_map.md`: 内容そのものは変更しないが、他 3 文書との
   記述(対応コマンド数・8 種ドーパント等)に矛盾が無いか照合。
5. `docs/tasks/README.md`: 作業終盤に `git log`/`git status` を再確認し、
   並行実行中の C-2 エージェントが同ファイルを更新していれば
   `git pull --rebase` の上で両者の変更をマージし、状態表が矛盾しないか
   通読確認(自分は C-2 の担当行(パリティ表の C-2 行、クラスタ校正関連)を
   変更しない)。
6. `docs/*.md` 全体を `未実装|未対応|未着手|TODO` で grep し、各ヒットが
   現在も正しい主張かを実装/テストと突き合わせて確認(誤りがあれば注記追加)。

## 4. テスト仕様

文書のみの変更のため自動テストは無いが、事実確認として以下を実行し記録する:

- `cmake --build build -j$(nproc)` が成功すること(ビルドが壊れていないことの前提確認)。
- `ctest --test-dir build -R sprocess_parity -V` の出力サマリ(`N/3 parity checks
  passing` 相当行)で、README/gap_analysis に書く「未達項目数」が実測と一致することを確認。
- `ctest --test-dir build -R "golden_flows|test_implant_materials|test_dual_pearson|test_oxidation"`
  が PASS することを確認し、W-7/W-8/C-1/C-3 の「解消済み」記載の裏付けとする。
- `git log --oneline | grep -i "C-2:"` をコミット直前に再実行し、C-2 の記載を
  最新化。

## 5. 完了条件 (DoD)

- [ ] `sentaurus_gap_analysis_v2.md` の各表・§4 が「解消済み/残存」を一目で判別できる。
- [ ] `structure_model_root_cause.md` §2 に是正済み注記あり。
- [ ] `IMPLEMENTATION_PLAN_v2.md` に W-1/W-3/W-4/W-7/W-8/C-1/C-3(+可能なら C-2)の
      完了状態表あり。
- [ ] `docs/sprocess_command_map.md` は無変更、他文書との矛盾なし。
- [ ] `docs/tasks/README.md` が C-2 エージェントの変更と矛盾なくマージされている
      (該当する場合)。
- [ ] grep で見つかる「未実装」等の記述が全て現状と一致(古い誤記ゼロ)。
- [ ] ビルド・該当テストが実際に PASS することを確認した上で記述している。

## 6. やらないこと

- `docs/tasks/*.md` 個別タスク仕様書の書き換え(歴史的記録として保持)。
- `tests/test_sprocess_parity.cpp` およびその CMake 登録の変更(C-2 の領域)。
- `docs/tasks/README.md` の C-2 パリティ行・クラスタ校正関連の記述変更。
- `src/` `tests/`(本タスク由来のもの以外) の変更。
- `sprocess_command_map.md` の再改訂(W-1 で完了済みのため)。
