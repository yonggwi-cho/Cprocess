# PA-5: MPI 分散 (長期・第一段階: MC 注入のみ)

## 目的

クラスタ実行に向けた MPI 対応の第一段階として、**MC 注入のみ**を
ランク分散する。MC はイオンごとに独立 (embarrassingly parallel) で
メッシュ分割が不要なため、最小の侵襲で線形スケーリングが得られる。
拡散ソルバーの領域分割 (METIS + ハロー交換) は明示的にスコープ外とし、
前提条件のみ本仕様の末尾に記す。

## 現状コード

- `src/mc_implant.cpp` — `apply_mc_implant`: イオンを `kChunk` 単位の
  チャンクに分け、チャンク k のシードは
  `splitmix64(p.seed ^ (0x9E3779B97F4A7C15ull*(k+1)))` で**チャンク番号
  のみから決まる** (スレッド数非依存・ビット再現)。この性質により、
  「ランク r がチャンク部分集合を担当し、ヒット配列を総和する」だけで
  1 ランク実行と**ビット一致**する分散が可能 (アモルファスモード)
- チャネリング + 損傷フィードバックはスレッド数/実行順で非決定的
  (コードコメントに明記済み) — MPI 分散では**損傷配列がランク間で
  共有されない**ため、チャネリング時は結果がランク数依存になる。
  仕様: チャネリング有効時は分散を許可するが、ドキュメントと
  ログに非決定性を明記する
- `CMakeLists.txt` — MPI 設定なし

## 実装手順

1. CMake:
   ```cmake
   option(CPROCESS_MPI "Enable MPI-distributed MC implant" OFF)
   if(CPROCESS_MPI)
     find_package(MPI REQUIRED COMPONENTS CXX)
     target_link_libraries(cprocess_core PUBLIC MPI::MPI_CXX)
     target_compile_definitions(cprocess_core PUBLIC CPROCESS_MPI)
   endif()
   ```
2. `src/mc_implant.cpp` の `apply_mc_implant` 末尾処理を拡張
   (`#ifdef CPROCESS_MPI` ガード):
   - 冒頭で `int rank=0, nprocs=1;` を取得
     (`MPI_Initialized` を確認し、未初期化なら rank=0/nprocs=1 の
     非分散動作 — ライブラリ側で MPI_Init は**しない**。呼び出し側
     アプリ/mpirun の責務とし、`cprocess` CLI の main に
     `#ifdef CPROCESS_MPI` の Init/Finalize を追加する)
   - チャンクループを `for (k = rank; k < nchunks; k += nprocs)` の
     ラウンドロビン割当に変更 (OpenMP 並列はランク内でそのまま有効。
     チャンクシードは k のみに依存するため割当方式は結果に影響しない)
   - 集約: ループ後に
     ```cpp
     MPI_Allreduce(MPI_IN_PLACE, hits_total.data(), nc,
                   MPI_UINT64_T, MPI_SUM, MPI_COMM_WORLD);
     ```
     同様に ChunkStat の 6 カウンタ + dsum/d2sum を 1 本の
     double/int64 バッファに詰めて Allreduce (順序独立な総和のみ。
     dsum/d2sum の加算順が変わるため Rp/dRp は 1 ランク実行と
     最終桁が異なりうる — ヒット配列 (整数) は**厳密一致**する。
     この差を仕様として明記)
   - 損傷配列 `dmg` も Allreduce (チャネリング無効時のみ意味が厳密)
3. デッキ/Python は無変更 (mpirun で外部起動)。pybind ビルドは
   CPROCESS_MPI=OFF を推奨とドキュメント (README に 1 行追記)
4. テスト: `tests/test_mpi_mc.cpp` 新設。CMakeLists では
   `if(CPROCESS_MPI)` 内で登録し、
   `add_test(NAME mpi_mc COMMAND ${MPIEXEC_EXECUTABLE} -n 4 test_mpi_mc)`

## テスト仕様

CPROCESS_MPI=ON + 4 ランク (`mpirun -n 4`):

1. **ビット一致 (アモルファス)**: 4×4×20 箱メッシュ、B 50keV 4e4 イオン
   channeling=false。4 ランクの Allreduce 後ヒット配列 (rank0 で検証) が、
   同一パラメータの 1 ランク実行 (テスト内で nprocs=1 相当の
   逐次チャンクループを直接回して再計算) と**全セルで整数一致**
2. **統計一致**: deposited/backscattered/transmitted の合計が 1 ランク値と
   一致 (整数)、Rp の差 < 1e-12 相対
3. **濃度場**: 変換後の conc 配列が 1 ランク再計算と相対差 < 1e-12

CPROCESS_MPI=OFF (既定・CI):
4. 全既存テストが無変更で PASS (ガード検証)。`MPI_Initialized`
   フォールバックにより OFF ビルドはコードパス完全同一

## 完了条件 (DoD)

- [ ] CPROCESS_MPI=OFF で全既存テスト PASS・結果無変更
- [ ] 4 ランクでヒット配列が 1 ランクと整数一致 (チャネリング無効時)
- [ ] ライブラリは MPI_Init しない (Initialized チェックで安全)
- [ ] チャネリング時の非決定性がログ/ドキュメントに明記されている
- [ ] コミットメッセージに `PA-5` を含める

## やらないこと

- 拡散ソルバーの領域分割。前提条件 (将来仕様の入力): Mesh への
  ゴーストセル層、分散 CSR とハロー交換、ILU の代替 (ランク毎
  block-Jacobi + 局所 ILU)、METIS 依存の導入判断
- 損傷フィードバックのランク間同期 (チャネリング決定性)
- 動的負荷分散 (ラウンドロビンで十分)
- Python からの MPI 起動サポート
