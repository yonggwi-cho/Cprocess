# P3-g: 1D/2D 高速モード

## 目的

較正・パラメータスイープ用途では横方向構造のない 1D (nx=ny=1) /
2D (ny=1) 実行で十分であり、3D フルメッシュ比 50〜100 倍の高速化が
見込める。本タスクは (1) `make_box_mesh` の nx=1/ny=1 動作を検証・保証、
(2) MC 注入の横方向脱出をオプトインの周期境界 (`lateral_wrap`、既定 off)
で折り返し可能に、(3) `Simulation.mesh1d(z, nz)` 便宜メソッドを追加する。
拡散ソルバーは FVM で次元非依存のため**変更なし** (P3_overview.md の
確定事項)。

前提タスク: なし (独立)。

## 現状コード

- `src/box_mesh.cpp` L13–17 — `make_box_mesh` の制約は `nx,ny,nz >= 1` のみ
  で、**nx=1 / ny=1 は既に構文上許可されている**。Kuhn 6-tet 分割
  (L31–41) は各 hex を主対角 v0–v7 共有の 6 tet に切るだけで、分割の
  整合性は隣接 hex 間の共有 quad 対角にしか依存しない。したがって
  **1 セル幅 (nx=1) でも退化せず成立する** — 頂点数は x 方向 2 枚、
  hex は 1 列になるだけで、tet 体積・面分類 (L48–65; 1 セル幅でも
  各境界面の全ノードは単一平面上にあり patch 分類は正常) とも問題ない。
  コード読解の結論: **プリズム分割の追加は不要**。overview の
  「成立するか検証 — しない場合は 1D 専用プリズム分割」の分岐は
  「成立する」側で確定し、成立をテストでロックする (下記テスト 1)。
  注意点: nx=1 で横 extent を深さと同程度に取ると tet が細長くなるが、
  1D 用途では横 extent をセル高さ dz 程度に選べばアスペクト比 ~1 に
  保てる (mesh1d の既定値で対処、手順 4)。
- `src/mc_implant.cpp` L294–311 `WalkParams` (`wrap`, `sx1..sy2`)、
  L320–323 `wrap1` ラムダ、L336–342 — 横方向脱出処理は既に存在する:
  `w.wrap` が真なら x/y を周期折返し、偽なら `Fate::out_of_domain`。
  L439 `w.wrap = !p.has_window;` — つまり**全面ビームでは既に周期
  wrap、window モードでのみ out_of_domain になる**。P3-g の
  `lateral_wrap` は「window モードでも wrap を強制する」オプトインとして
  実装する (既定 off = 現行挙動)。なお wrap1 は fmod 折返し
  (トーラス) であり鏡映ではない — 一様横方向統計では等価。
- `include/cprocess/mc_implant.hpp` L51–72 `McImplantParams` —
  `has_window`/`x1..y2` はあるが wrap 制御フラグはない。
- `include/cprocess/process.hpp` L56–62 / `src/process.cpp` L303–401 —
  `proc::implant_mc(st, species, dose, energy_kev, ions, tilt, rotation,
  seed, threads, channeling, has_window, x1, x2, y1, y2, seed_damage, log)`。
- `python/_cprocess.cpp` L443–461 — `proc_implant_mc` バインディング
  (`py::arg("damage") = false` まで)。L380–392 — `proc_mesh_box`。
- `python/cprocess/simulation.py` L92–99 — `Simulation.mesh(x, y, z, nx, ny,
  nz, x0=..., y0=..., z0=...)` が `proc_mesh_box` を µm→cm 変換して呼ぶ。
  L124–160 — `Simulation.implant(..., window=None, ...)` が
  `proc_implant_mc` を呼ぶ。1D 便宜メソッドは存在しない。
- `python/test_comprehensive.py` — Python 層の回帰テスト。

## 実装手順

1. **メッシュ**: `make_box_mesh` 本体は無変更 (nx=1 は既に動く)。
   保証はテストで行う (テスト 1)。ドキュメントコメント (box_mesh.cpp
   L9–12) に「nx=1/ny=1 (1D/2D モード) をサポートする」旨を 1 行追記。
2. **MC lateral_wrap**: `McImplantParams` (mc_implant.hpp L51–72) に
   ```cpp
   bool lateral_wrap = false;  // wrap lateral exits even in window mode
   ```
   を追加し、`src/mc_implant.cpp` L439 を
   ```cpp
   w.wrap = !p.has_window || p.lateral_wrap;
   ```
   に変更。`walk_ion` は無変更 (wrap 機構は既存)。既定 false で全経路
   現行とビット一致 (RNG 消費・分岐とも不変)。
3. **proc:: 層** (CLAUDE.md 必須ラッパー則の対象):
   `proc::implant_mc` (process.hpp L56–62, process.cpp L303) の末尾
   `seed_damage` の直前に `bool lateral_wrap = false` 引数を追加し、
   `p.lateral_wrap = lateral_wrap;` を設定 (process.cpp L326 付近)。
   physical-resist 分岐 (L336 `p.has_window = false;`) では wrap は
   もともと有効なので追加処理不要。
4. **pybind** (`python/_cprocess.cpp`):
   - `proc_implant_mc` (L443) のラムダと `py::arg` 列に
     `py::arg("lateral_wrap") = false` を `damage` の前に追加
     (C++ 引数順と一致させる)。
   - `mesh1d` は新規 C++ 関数を追加しない (既存 `proc::mesh_box` の
     退化呼び出しで実現するため、新規バインディングは
     `lateral_wrap` の 1 個のみ)。
5. **Simulation メソッド** (`python/cprocess/simulation.py`):
   - `implant()` (L124) のキーワードに `lateral_wrap: bool = False` を
     追加し、mc パスで `proc_implant_mc(..., bool(lateral_wrap),
     bool(damage))` と渡す (analytic パスでは無視、mc=False かつ
     lateral_wrap=True なら ValueError)。
   - `mesh()` の直後に便宜メソッドを追加:
     ```python
     def mesh1d(self, z: float, nz: int,
                lateral: float = 0.0) -> "Simulation":
         """1D fast-mode mesh: nx=ny=1, depth z [um], nz cells.

         `lateral` is the x/y extent [um]; default 0 -> one cell height
         z/nz (keeps tet aspect ratio ~1). Delegates to mesh().
         """
         lat = lateral if lateral > 0 else z / nz
         return self.mesh(x=lat, y=lat, z=z, nx=1, ny=1, nz=nz)
     ```
     `self.mesh()` 経由なので単位変換・`_emit`・`return self`
     (メソッドチェーン) は既存実装をそのまま流用する。
6. 拡散 (`diffusion.cpp`)・deposit/etch・oxidize は**変更しない**。
   1D メッシュでもセル/面ループは次元非依存に動く。
7. 調整可能定数は本タスクでは発生しない。将来 lateral 既定比率等を
   可変にする場合は P1-10 の `ParamDB`
   (`include/cprocess/param_db.hpp`) キー (例 `"mesh.lat_aspect"`) で
   読むこと。ハードコードの新定数は追加しない。

## テスト仕様

新規 `tests/test_fast_1d2d.cpp` (main + `test_util.hpp`)。
`CMakeLists.txt` L76 の `foreach(t ...)` に `fast_1d2d` を追加。

1. **nx=1 メッシュ健全性**: `make_box_mesh(0,0.01e-4, 0,0.01e-4, 0,0.5e-4,
   1,1,50)` について (a) `cells.size() == 6*50`、(b) 全 `cell_vol[i] > 0`
   かつ `Σ cell_vol == 0.01e-4 * 0.01e-4 * 0.5e-4` (相対誤差 < 1e-12)、
   (c) 全境界面の patch が 0..5 に分類 (例外なく finalize 完了)、
   (d) ny=1 のみの 2D ケース `(8,1,50)` も同様に (a)–(c)。
2. **1D vs 3D 深さプロファイル < 3%**: 解析 Gaussian 注入
   (`apply_implant`, B, Rp=50nm, dRp=20nm, dose 1e14) + 拡散
   (1000°C=1273.15K, 600s) を
   (a) 1D: `make_box_mesh(..., 1, 1, 50)`、
   (b) 3D: 同 z 範囲で 8×8 lateral (`8, 8, 50`)
   の両方で実行。z を 50 層にビン分けし、各層で体積加重平均濃度を計算
   ((b) は x,y 平均)。ピーク濃度の 1e-3 倍以上の層について
   相対差 < 3% を assert。
3. **≥50 倍高速化**: テスト 2 と同一物理の拡散ソルブ (実装は
   `DiffusionSolver` 直呼びで可) を threads=1
   (`omp_set_num_threads(1)`) で `std::chrono::steady_clock` 計測し、
   `t_3d >= 50.0 * t_1d` を assert。セル数比は 64 倍でソルバーは
   超線形なので余裕があるが、ノイズ対策に 1D 側は 5 回計測の最小値、
   拡散時間は 3D 側が数秒オーダーになるよう nz・時間を調整する
   (nz=100, 10 ステップ以上を目安)。
4. **lateral_wrap dose 損失 < 0.1%**: 狭い 1D 風メッシュ
   (0.02×0.02×0.5 µm, nx=ny=2) で MC 注入 (B 30keV, ions=100000,
   channeling=false, seed=3, window=(全域)) を
   (a) `lateral_wrap=false`: `stats.out_of_domain > 0` を確認
   (現行挙動: window モードで横脱出)、
   (b) `lateral_wrap=true`: `stats.out_of_domain == 0` (受入基準
   < 0.1% を厳密化) かつ
   `deposited+backscattered+transmitted+in_mask+unbinned == ions`。
   また has_window=false では lateral_wrap の真偽で結果がビット一致
   (`w.wrap` が両方 true になるため) であることも assert。
5. **既存回帰**: `test_mesh` / `test_mc` / `test_diffusion` /
   `test_integration` / `test_resist` を含む既存全テストが無変更で PASS
   (lateral_wrap 既定 off は全経路ビット不変)。

Python (`python/test_comprehensive.py` に追加、CLAUDE.md 必須):

- `Simulation().mesh1d(z=0.5, nz=50)` → `sim._st.mesh` のセル数 300、
  メソッドチェーンで `.init(...).implant("B", 1e14, energy=30)` が通る。
- `mesh1d` + analytic implant + `diffuse` の depth プロファイル取得が
  例外なく動く (値は C++ テスト 2 でカバー済みのため存在確認レベル)。
- `implant("B", 1e14, energy=30, mc=True, ions=20000,
  window=(0, 0.01, 0, 0.01), lateral_wrap=True)` の返す
  `ImplantResult.stats.out_of_domain == 0`。
- `mc=False, lateral_wrap=True` が ValueError。

## 完了条件 (DoD)

- [ ] nx=1/ny=1 の box mesh がテストで保証 (プリズム分割は追加しない —
      Kuhn 分割が 1 セル幅で成立することをコード読解で確認済み)
- [ ] `McImplantParams::lateral_wrap` (既定 false) と
      `w.wrap = !has_window || lateral_wrap` の 1 行変更
- [ ] `proc::implant_mc` に `lateral_wrap` 引数、pybind
      `proc_implant_mc` に `py::arg("lateral_wrap")=false`、
      `Simulation.implant(lateral_wrap=)` と `Simulation.mesh1d()` 追加
- [ ] `tests/test_fast_1d2d.cpp` の 5 テスト PASS、CMakeLists 登録済み
- [ ] `python/test_comprehensive.py` に mesh1d / lateral_wrap テスト追加、
      PASS (`build/_cprocess*.so` を `python/cprocess/` へコピーして実行)
- [ ] 既存全テスト PASS (`ctest --test-dir build`)、lateral_wrap 既定 off
      で全結果ビット不変
- [ ] コミットメッセージに `P3-g` を含める

## やらないこと

- 拡散ソルバー・酸化・応力の 1D 特化 (FVM は次元非依存 — 変更なし、
  overview の確定事項)
- 1D 専用プリズム分割 (Kuhn 分割が nx=1 で成立するため不要)
- 鏡映 (specular) 境界 — wrap は既存の fmod 周期折返しを流用
- lateral_wrap の text deck (deck.cpp) コマンド化 (Python 層が公式
  入口; 必要なら別タスク)
- 2D 専用の可視化・出力形式 (VTU は退化 3D として出力される)
- チャネリング on の MC 決定性改善 (現行仕様のまま)
