# M-5: 品質駆動リメッシュ (ALE 後の自動修復パス)

## 目的

`ale_move` (src/ale_mover.cpp) は界面ノードを法線方向へ動かすが、変位が
セル寸法に近づくと sliver (q → 0) や反転寸前のセルが残る。反転ガードは
「動かさない」だけなので表面が乱れる。酸化統合 (P1-6) では界面が毎ステップ
動くため、**移動後にメッシュ品質を自動回復するパス**が必須。
平滑化 → 分割 の 2 段修復を実装する。

## 現状コード

- `src/remesh.cpp` — `tet_quality` (正規 tet=1 の品質)、
  `mesh_quality(m, sliver_thresh=0.1)` (min_q / mean_q / n_sliver / worst_cell)、
  `laplacian_smooth(m, iters, omega)` (内部ノードのみ、反転ガード付き)
- M-1 完了後: `split_edges` (finalize 内包、cell_parent 付き) が使える
- `src/ale_mover.cpp` — `ale_move(...)` は境界/界面ノードを動かし
  `AleResult{n_moved, n_skipped}` を返す。finalize は呼び出し側
- `include/cprocess/field_transfer.hpp` — `transfer_field_nearest` (使わない。
  M-1 の `redistribute_field` で足りる)

## 実装手順

1. `include/cprocess/remesh.hpp` に追加:
   ```cpp
   struct RepairResult {
     int n_smoothed = 0;     // 平滑化で動いたノード数 (全パス合計)
     int n_split = 0;        // 分割されたエッジ数
     double min_q_before = 0, min_q_after = 0;
     std::vector<int> cell_parent;  // 分割があった場合の親マップ
                                    // (分割なしなら iota)
   };
   // 品質 q < q_thresh のセルを修復する。手順:
   //   1) laplacian_smooth を smooth_iters 回
   //   2) まだ q < q_thresh のセルについて、その最長エッジを split_edges で分割
   //   3) もう一度 laplacian_smooth
   // を max_rounds 回まで繰り返し、n_sliver == 0 になったら早期終了。
   RepairResult repair_quality(Mesh& m, double q_thresh = 0.1,
                               int max_rounds = 3, int smooth_iters = 5);
   ```
2. 実装 (src/remesh.cpp):
   - (a) `min_q_before = mesh_quality(m, q_thresh).min_q` を記録
   - (b) 各 round:
     - `laplacian_smooth(m, smooth_iters, 0.5)` → smooth 後に
       `m.finalize()` (cell_vol/cell_cent がノード移動で古くなるため。
       laplacian_smooth 自体は nodes しか変更しない)
     - `mesh_quality` を再計算。`n_sliver == 0` なら break
     - sliver セル (q < q_thresh) を列挙し、各セルの**最長エッジ**
       (kEdges 6 通りから) を候補に集める。重複エッジは除去
       (両端ノード番号を正規化した pair の set)
     - `split_edges(m, candidates)` を呼ぶ。`cell_parent` を
       ラウンドをまたいで**合成**する: `parent_total[i] =
       parent_total_prev[result.cell_parent[i]]`
   - (c) `min_q_after` を記録して返す
   - 注意: 分割は sliver を必ず直すとは限らない (针状 tet)。max_rounds で
     打ち切り、結果を戻り値で報告するだけにする (throw しない)
3. `ale_move` との結線はしない (P1-6 のスコープ)。ここでは独立関数として
   実装し、テストで「ALE → repair」の連続適用を検証する

## テスト仕様

`tests/test_remesh.cpp` に追加:

1. **人工 sliver の修復**: 6×6×6 箱メッシュの内部ノード 1 個を、隣接セルが
   q < 0.05 になるまで意図的にずらす (座標を辺長の 0.45 倍だけ移動し
   finalize)。`repair_quality(m, 0.1)` 後に
   `mesh_quality(m).min_q > 0.05` かつ `min_q_after > min_q_before`
2. **健全メッシュは無変更**: 生成直後の箱メッシュ (min_q ≈ 0.3+) に対し
   `repair_quality` → `n_split == 0` かつノード座標が不変 (bit 一致)
   ※ laplacian_smooth は品質が q_thresh 以上なら走らせない実装でもよいし、
   round 冒頭で `n_sliver == 0` → 即 return でもよい (後者を推奨)
3. **ALE 連続適用**: 4×4×8 箱メッシュで zmax 界面を
   `ale_move` (v_n·dt = セル高の 0.3 倍) → `finalize` → `repair_quality`
   を 5 回繰り返す。各回で `min_q > 0.02` を維持し、最終メッシュで
   体積が解析値 (元の体積 + 移動量×断面積) と相対差 < 5% で一致
4. **cell_parent 合成**: テスト 1 で fields を模した一様配列を
   `redistribute_field(conc, result.cell_parent)` で転写 → 質量保存 < 1e-12

## 完了条件 (DoD)

- [ ] `repair_quality` 実装 (throw しない、結果は戻り値報告)
- [ ] 上記 4 テスト追加、全テスト PASS
- [ ] ラウンドをまたいだ `cell_parent` の合成が正しい (テスト 4)
- [ ] コミットメッセージに `M-5` を含める

## やらないこと

- エッジ縮約・2-3/3-2 フリップ (M-3 / M-4)。分割+平滑化で直らない sliver は
  報告のみ
- `ale_move` 内部への結線 (P1-6 で酸化ステップが repair を呼ぶ)
- proc:: / Python API (P1-6 の oxidize が内部で使う)
- 並列化 (PA-1)
