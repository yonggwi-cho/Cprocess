# M-2: 適応細分化 (勾配駆動リファイン、P1-8 と同一タスク)

## 目的

急峻な濃度プロファイル (TED 初期の implant テール、界面近傍) を粗い初期メッシュ
のまま解くと数値拡散で崩れる。M-1 の `split_edges` を基盤に、**濃度勾配が大きい
面のセル辺を自動選択して分割する適応細分化ループ** `proc::refine` を実装し、
デッキ/Python から 1 ステップとして呼べるようにする。IMPLEMENTATION_PLAN の
P1-8 (適応リメッシュ最小限) と M-2 は同一タスク。

## 現状コード

- `src/remesh.cpp` / `include/cprocess/remesh.hpp` —
  `split_edges(Mesh&, edges)`: 競合スキップ付き一括分割、`cell_parent` を返し
  内部で `finalize()` 済み。`redistribute_field(conc, cell_parent)`:
  intensive 場の親→子コピー転写。`repair_quality(Mesh&, fields*, ...)`:
  fields ポインタ配列を受け取り分割/フリップを通して場を転写する overload あり
- `src/process.cpp` — `SimState::fields` は
  `std::map<std::string, std::vector<double>>` (セル濃度 cm^-3)。
  `need_mesh(st)` / `fmt(...)` などのヘルパは無名 namespace 内
- `src/diffusion.cpp` — `DiffusionSolver::gradients` は面近傍の加重最小二乗
  勾配 (参考実装。本タスクでは**使わない** — 下記の面ベース指標で十分)
- `Mesh::faces` — 各 `Face` は `owner`/`neigh` (neigh<0 は境界)、
  `Mesh::cell_vol`/`cell_cent` は finalize() 後に有効
- 分割エッジの自動選択ロジックは存在しない (M-1 の「やらないこと」)

## 実装手順

1. `include/cprocess/remesh.hpp` に内部 API を追加:
   ```cpp
   struct RefineResult {
     int n_passes = 0;        // 実行したパス数
     int n_split_total = 0;   // 分割エッジ総数
     int n_cells_before = 0, n_cells_after = 0;
   };
   // 勾配指標に基づく適応細分化。conc はセル濃度 (呼び出し時点のセル数)。
   // fields に列挙された全場 (conc 自身を含む) を各パスの cell_parent で
   // redistribute_field 転写する。max_growth: セル数が
   // n_cells_before * max_growth を超えたら以降のパスを打ち切る。
   RefineResult refine_gradient(Mesh& m, const std::vector<double>& conc_key,
                                std::vector<std::vector<double>*>* fields,
                                double rel_grad_thresh, int max_passes,
                                double max_growth = 4.0);
   ```
   注意: `conc_key` は指標計算のキー種のみ。実装内では fields に含まれる
   キー種のポインタを毎パス参照する (下記 2(f) 参照)。
2. `refine_gradient` の実装 (src/remesh.cpp)。各パス:
   - (a) `global_max = max_i conc[i]` を計算 (conc はキー種の現在値)。
     `global_max <= 0` なら候補ゼロとして終了
   - (b) **面を走査** (エッジではない): 各内部面 `f` (`f.neigh >= 0`) について
     指標を評価する。owner/neigh のセル値を `C_o`, `C_n` として、
     ```
     |C_o - C_n| > rel_grad_thresh * max(C_o, C_n, 1e-3 * global_max)
     ```
     を満たす面を「急勾配面」とする (床値 1e-3·global_max はゼロ近傍の
     比率発散を防ぐ)
   - (c) 急勾配面ごとに、**owner セルの最長エッジ** (6 辺のうちノード間距離
     最大のもの、`kEdges` テーブル流用) を分割候補に加える。ノード対は
     (a<b) に正規化し `std::set` で重複排除
   - (d) 候補が空ならループ終了。空でなければ `split_edges(m, candidates)`
     を呼ぶ (競合スキップは split_edges 側の仕様に任せる。スキップ分は
     次パスで再評価される)
   - (e) fields の全ベクトルを `redistribute_field(*v, sr.cell_parent)` で
     転写する (conc_key が fields 内のベクトルを指していればここで一緒に
     更新される — 呼び出し側はそう渡すこと)
   - (f) 次パスの指標評価には転写後のキー種を使う。実装上は
     `refine_gradient` の第 2 引数を「fields 内のキー種のインデックス」
     `int key_index` に変更してもよい — **推奨シグネチャ**:
     ```cpp
     RefineResult refine_gradient(Mesh& m,
                                  std::vector<std::vector<double>*>& fields,
                                  int key_index, double rel_grad_thresh,
                                  int max_passes, double max_growth = 4.0);
     ```
     (以降この形を正とする。key_index が fields 範囲外なら
     `std::invalid_argument`)
   - (g) セル数が `n_cells_before * max_growth` 以上になったら打ち切り。
     `max_passes` パスで終了
3. `include/cprocess/process.hpp` / `src/process.cpp` に proc:: 関数を追加:
   ```cpp
   // species の濃度勾配が急峻な領域を適応細分化する。
   // rel_grad_thresh: 面をまたぐ相対濃度差の閾値 (無次元, 既定 0.5)
   // max_passes: 分割パス数上限 (既定 2)
   void refine(SimState& st, const std::string& species,
               double rel_grad_thresh, int max_passes,
               std::ostream* log = nullptr);
   ```
   実装:
   - `need_mesh(st)`; `st.fields` に species が無い/空なら
     `std::runtime_error("refine: no field '<species>'")`
   - `st.fields` の**全**ベクトルへのポインタを `std::vector<std::vector<double>*>`
     に集め (map の走査順で)、species のインデックスを key_index として
     `refine_gradient(st.mesh, fields, key_index, rel_grad_thresh, max_passes)`
     を呼ぶ (max_growth は既定 4.0 のまま)
   - `st.has_stack` が true なら stack は分割対象外 (触らない)。
     ログに `[refine] <species> passes=N split=M cells X -> Y` を出力
4. pybind11 バインディング (python/_cprocess.cpp):
   ```cpp
   m.def("proc_refine",
       [](SimState& st, const std::string& species, double thresh, int passes) {
         std::ostringstream log;
         proc::refine(st, species, thresh, passes, &log);
         return log.str();
       },
       py::arg("state"), py::arg("species"),
       py::arg("threshold") = 0.5, py::arg("passes") = 2);
   ```
5. `Simulation` メソッド (python/cprocess/simulation.py)。単位変換なし
   (閾値は無次元):
   ```python
   def refine(self, species: str, threshold: float = 0.5,
              passes: int = 2) -> "Simulation":
       """Adaptively split mesh edges where `species` has steep gradients."""
       self._emit(_c.proc_refine(self._st, species, float(threshold),
                                 int(passes)))
       return self
   ```
6. テスト追加: `tests/test_remesh.cpp` に C++ ケース、
   `python/test_comprehensive.py` に Python ケース。CMakeLists.txt の
   foreach は既存 `remesh` を使うため変更不要

## テスト仕様

`tests/test_remesh.cpp` に追加 (SimState を使うテストは
`tests/test_integration.cpp` でもよいが、以下は refine_gradient を直接叩く
形と proc::refine を叩く形を 1 本ずつ):

1. **高勾配帯の細分化**: 箱メッシュ 0.4×0.4×0.8 µm (cm 単位)、
   nx=ny=6, nz=16。B 場に縦 Gaussian
   `C(z) = 1e19 * exp(-(z_top - z - Rp)^2 / (2 dRp^2))`
   (Rp=0.1µm=1e-5 cm, dRp=0.02µm=2e-6 cm) をセル centroid で評価して設定。
   `proc::refine(st, "B", 0.5, 2)` 実行後:
   - 深さ帯 Rp±2dRp に centroid を持つセル数が実行前の **> 1.5 倍**
   - 遠方帯 (z_top − z > Rp + 6dRp) のセル数は実行前と**同数** (不変)
2. **ドーズ保存**: テスト 1 の前後で Σ C_B·V の相対差 **< 1e-9**
   (redistribute_field は厳密コピーなので実際は < 1e-12 のはず)。
   B 以外にもう 1 種 (P 一様 1e15) を入れておき、P の質量も < 1e-9 で保存
3. **成長上限**: rel_grad_thresh=0.01, max_passes=10 で実行しても
   `n_cells_after <= 4 * n_cells_before` (max_growth=4 が効く)
4. **一様場は無変化**: 一様 1e15 の場で refine → `n_split_total == 0`、
   セル数不変
5. **品質**: テスト 1 実行後 `mesh_quality(m).min_q > 0` (反転なし)、
   `total_volume()` 相対差 < 1e-12
6. 既存の remesh / ale / integration テストがすべて PASS

`python/test_comprehensive.py` に追加:

- `sim.mesh(0.4, 0.4, 0.8, 6, 6, 16)` → `sim.implant("B", 1e14, energy=30)`
  → `n0 = sim.n_cells` → `sim.refine("B", threshold=0.5, passes=2)` →
  `sim.n_cells > n0` かつ `sim.dose("B")` が refine 前後で相対差 < 1e-9。
  メソッドチェーン (`sim.refine(...)` が sim を返す) を確認

## 完了条件 (DoD)

- [ ] `refine_gradient` (remesh.hpp/cpp) と `proc::refine` (process.hpp/cpp) 実装
- [ ] `proc_refine` バインディング + `Simulation.refine` + Python テスト追加
- [ ] 上記 C++ 5 テスト + Python 1 テスト追加、
      `ctest --test-dir build` 全 PASS、`python3 python/test_comprehensive.py` PASS
- [ ] コミットメッセージに `M-2` を含める

## やらないこと

- ノード補間勾配・最小二乗勾配による指標 (面 2 セル差で十分。
  `DiffusionSolver::gradients` は流用しない)
- 粗視化 (M-3)、異方性指標 (M-6)、界面距離ベースの指標
- 拡散ソルバー内からの自動トリガ (ユーザが工程間で明示的に refine を呼ぶ)
- stack メッシュ (photo 中) の細分化
- 並列化 (PA-1)
