# PA-4: GPU オフロード (長期・第一段階)

## 目的

10⁷ セル級の 3D フル CMOS シミュレーションに向け、反復ソルバーの
ホットカーネル (SpMV / AXPY / 内積) を GPU にオフロードする。
第一段階は **OpenMP target offload** を採用する (CUDA/HIP ではなく)。
理由: 単一ツールチェーン (gcc/clang) で CPU フォールバックと共存でき、
既存の `#pragma omp` コードベースと整合するため。性能を極めるのは
後続タスク (CUDA 専用カーネル化) とし、本タスクは「正しく動く
オフロード経路と CPU 完全互換のフォールバック」を確立する。

## 現状コード

- `src/sparse.cpp` — `CSR::mul`, `dotv`, `axpy`, `inv_diag` は
  `#pragma omp parallel for` (ホスト並列) 済み。cg/bicgstab は
  これらの組合せ
- `ILU0::apply` のレベル三角ソルブはレベル内並列だが逐次依存が強く
  GPU 向きでない — 本タスクでは **GPU 経路は Jacobi 前処理のみ**とする
- `CMakeLists.txt` — オフロード関連の設定なし

## 実装手順

1. CMake オプション追加:
   ```cmake
   option(CPROCESS_GPU "Enable OpenMP target offload" OFF)
   if(CPROCESS_GPU)
     target_compile_definitions(cprocess_core PUBLIC CPROCESS_GPU)
     # gcc: -foffload=nvptx-none / clang: -fopenmp-targets=nvptx64
     # 検出はせず、キャッシュ変数 CPROCESS_OFFLOAD_FLAGS をユーザ指定とする
     target_compile_options(cprocess_core PUBLIC ${CPROCESS_OFFLOAD_FLAGS})
     target_link_options(cprocess_core PUBLIC ${CPROCESS_OFFLOAD_FLAGS})
   endif()
   ```
2. `src/sparse.cpp` に GPU 版 CG を追加 (`#ifdef CPROCESS_GPU` ガード):
   ```cpp
   SolveResult cg_jacobi_gpu(const CSR& A, const std::vector<double>& b,
                             std::vector<double>& x, double rtol, int maxit);
   ```
   - ソルバー呼び出しスコープでデバイスバッファを 1 回だけ確立:
     ```cpp
     #pragma omp target data map(to: ptr[0:n+1], col[0:nnz], val[0:nnz], \
                                     b_[0:n], M[0:n])                    \
                             map(tofrom: x_[0:n])                        \
                             map(alloc: r[0:n], z[0:n], p[0:n], q[0:n])
     ```
     (std::vector は `.data()` の生ポインタで map する。スカラー
     リダクション (dot) は `map(tofrom:)` したスカラーに
     `#pragma omp target teams distribute parallel for reduction(+:s)`)
   - 各カーネル: `#pragma omp target teams distribute parallel for`
     で SpMV / AXPY / Jacobi 適用を実装。反復ループ自体はホスト側
     (収束判定のためドット積結果のみホストへ戻す)
   - **CPROCESS_GPU 未定義時はこの関数を宣言しない** (ヘッダも
     `#ifdef` ガード)。既存の `cg_jacobi` 等は一切変更しない
3. `DiffusionSolver::run` の切替は行わない (本タスクはソルバー基盤のみ。
   拡散への結線は GPU 実機での性能検証後の後続タスク)
4. テスト: `tests/test_gpu.cpp` を新設し、**CMakeLists では
   `if(CPROCESS_GPU)` の中でのみ** `foreach` と同形式で登録する

## テスト仕様

`tests/test_gpu.cpp` (CPROCESS_GPU=ON ビルドでのみコンパイル/実行):

1. **正しさ**: 2D 5 点ラプラシアン (50×50, test_solver の `laplacian2d` を
   共有ヘッダ化するか複製) で `cg_jacobi_gpu` と `cg_jacobi` の解の
   最大差 < 1e-9
2. **収束**: 反復回数が CPU 版の ±20% 以内
3. **ゼロ RHS**: x=0 を返し converged=true

CPROCESS_GPU=OFF (既定・CI 相当):
4. 全 15+ 既存テストが**無変更で** PASS すること (ガードの検証)。
   これが本タスクの主要な受入条件 — オフロード実機がない環境でも
   マージ可能であること

GPU 実機がない開発環境では、gcc のホストフォールバック
(`-foffload=disable` 相当でホスト実行) でテスト 1-3 を通す。
仕様: `omp_get_num_devices()==0` でもテストは正しさのみ検証するので
成立する (target リージョンはホストで実行される)。

## 完了条件 (DoD)

- [ ] CPROCESS_GPU=OFF で全既存テスト PASS (バイナリ・結果とも無変更)
- [ ] CPROCESS_GPU=ON がビルド可能 (少なくともホストフォールバックで
      test_gpu PASS)
- [ ] `cg_jacobi_gpu` の解が CPU と < 1e-9 で一致
- [ ] コミットメッセージに `PA-4` を含める

## やらないこと

- CUDA/HIP ネイティブカーネル (後続タスク)
- ILU 前処理の GPU 化 (Jacobi のみ)
- DiffusionSolver / MC への結線 (性能実測後の判断)
- マルチ GPU、ユニファイドメモリ最適化
