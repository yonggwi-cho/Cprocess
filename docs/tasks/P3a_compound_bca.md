# P3-a: 真の化合物 BCA (成分元素ごとの二体衝突)

## 目的

現行 MC 注入 (`src/mc_implant.cpp`) は SiO₂ 等の化合物を Bragg 則の
単一有効元素 (SiO₂: Z=10, M=20) で近似しており、軽元素 (O) との衝突で
大きく散乱される軽イオン (B, H) の後方散乱・飛程分布を再現できない。
本タスクでは `TargetMaterial` を成分元素リスト `{Z_i, M_i, x_i}` に拡張し、
各衝突で衝突相手元素を数密度比 x_i に比例した確率で抽選する
「真の化合物 BCA」に置き換える。電子的阻止能は Bragg 則 (成分加重和) を
維持し、自由行程は総数密度から取る。**単一元素材料 (Si 等) は成分 1 個の
退化ケースとして現行実装とビット一致**でなければならない
(P3_overview.md の確定事項)。

前提タスク: なし。変更は `mc_implant.hpp` / `mc_implant.cpp` (+ テスト) のみ。
Python 表面 (pybind/Simulation) のシグネチャ変更は不要 — 材料表は
`process.cpp` の `st.mat_table` 経由で内部的に流れるだけである。

## 現状コード

- `include/cprocess/mc_implant.hpp` L14–20 — `TargetMaterial { name, z, m, n,
  crystal_si }`: 有効単一元素。L23–31 のファクトリ
  `target_silicon()/target_photoresist()/target_oxide()/target_vacuum()` が
  組込み材料を返す (実装は `src/mc_implant.cpp` L172–188。
  `target_oxide()` は Z=10, M=20, N=6.62e22 の有効原子)。
- `src/mc_implant.cpp` L123–153 — `ScatterTable`: sin²(θ/2) を換算量
  (ε, b) の 160×160 log-log 格子に表化。**ZBL 万能遮蔽関数は Z に依存しない
  換算量のみの関数なので、このテーブル自体は (Z1,Z2) に依存しない万能
  テーブルである** (Z 依存は `screening_length_cm` / `reduced_energy_per_ev`
  経由で ε, b への換算係数に入る)。L153 `scatter_table()` が
  function-local static で 1 個だけ構築、L465 で worker 起動前に
  `scatter_table();` を呼び構築レースを回避している。
- L269–292 — `MatConstants` / `make_mat_constants(z1, m1, t)`: 材料ごとに
  `flight = N^(-1/3)`, `pmax`, `inv_a`, `eps_per_ev`, `tmax_fac`,
  `mass_ratio`, `els_fac` (Lindhard-Scharff k×N×flight) を前計算。
  **(Z1,Z2) ペア固有量はすべてここに集約されている。**
- L315–398 — `walk_ion`: 1 イオンの輸送。RNG 消費順は固定:
  (1) 初期位置 x, y で `rng.u()` × 2 (L316)、以降ステップごとに
  (2) 衝突径数 `rng.u()` (L372)、(3) 方位角 `rng.u()` (L392)。
  衝突相手は常に材料の有効元素 1 種 (`mc_.eps_per_ev` 等を直接使用、
  L374–375)。
- L513–537 — チャンク決定性: チャンク k ごとに
  `Rng rng(splitmix64(p.seed ^ (0x9E3779B97F4A7C15ull*(k+1))))` で独立
  ストリームを張り、スレッド数によらずビット一致 (channeling off 時)。
- ダメージ生成 (L379–387) は `crystal_si` 材料のみ、閾値 `kEdSi=15eV`。
- `src/process.cpp` L419 — `st.mat_table = {target_silicon(),
  target_photoresist(), target_vacuum()}` (photo() が構築)。
- テスト: `tests/test_mc.cpp` が散乱積分・Rp 回帰・決定性を検証済み。

## 実装手順

1. `include/cprocess/mc_implant.hpp` の `TargetMaterial` に成分リストを追加:
   ```cpp
   struct TargetComponent {
     int z = 14;        // atomic number
     double m = 28.086; // atomic mass [amu]
     double x = 1.0;    // number fraction (sum over components == 1)
   };
   struct TargetMaterial {
     const char* name = "Si";
     int z = 14;            // Bragg-rule effective Z (kept: logging/compat)
     double m = 28.086;     // effective M [amu]
     double n = 4.99e22;    // TOTAL atomic density [cm^-3] (all components)
     bool crystal_si = false;
     std::vector<TargetComponent> comp;  // empty => single element {z, m, 1.0}
   };
   ```
   `comp.empty()` は「有効単一元素 = 成分 1 個 {z,m,1.0}」と同義とし、
   既存の集約初期化 `{"Si", kZt, kMt, kNt, true}` を無変更で通す
   (完全後方互換)。
2. `target_oxide()` (mc_implant.cpp L180–184) に成分を設定:
   ```cpp
   TargetMaterial t{"SiO2", 10, 20.0, 6.62e22, false};
   t.comp = {{14, 28.086, 1.0/3.0}, {8, 15.999, 2.0/3.0}};
   return t;
   ```
   `target_photoresist()` / `target_vacuum()` は comp 空のまま
   (有効単一元素継続 — overview の対象は SiO₂/Si₃N₄ 型化合物)。
   Si₃N₄ 用ファクトリ `target_nitride()` も追加:
   密度 3.1 g/cm³, 7 原子 / 140.28 g/mol → N=9.32e22,
   comp = {{14,28.086,3/7.}, {7,14.007,4/7.}}, 有効 Z=(3·14+4·7)/7=10,
   M=140.28/7≈20.04。
3. `MatConstants` を成分別に分割 (mc_implant.cpp L270–292):
   ```cpp
   struct CompConstants {   // per (Z1, Z2_i) pair
     double inv_a, eps_per_ev, tmax_fac, mass_ratio;
     const ScatterTable* table;  // cached per-(Z1,Z2) table (step 4)
     double x_cum;               // cumulative number fraction for the draw
   };
   struct MatConstants {
     double flight, pmax, els_fac;  // material-level (total N / Bragg sum)
     bool crystal_si; std::vector<CrystalAxis> chan_axes;
     std::vector<CompConstants> comps;  // size 1 => legacy fast path
   };
   ```
   `make_mat_constants`:
   - `flight = pow(t.n, -1./3.)`, `pmax = flight/sqrt(M_PI)` は総数密度から
     (現行式そのまま — 確定事項「自由行程は総数密度から」)。
   - 成分ごとに `inv_a = 1/mc::screening_length_cm(z1, z_i)`,
     `eps_per_ev = mc::reduced_energy_per_ev(z1, m1, z_i, m_i)`,
     `tmax_fac = 4 m1 m_i/(m1+m_i)²`, `mass_ratio = m1/m_i` を計算。
   - 電子的阻止能 (Bragg 則):
     `els_fac = N * flight / sqrt(1000) * Σ_i x_i * kls(z1, m1, z_i)` と
     成分加重和にする (kls は現行 L287–288 の Lindhard-Scharff 式を
     z2→z_i にした関数に括り出す)。**comp が空 (成分 1 個) の場合は
     現行 L277–291 と全く同じ演算列・同じ丸めで計算すること**
     (加重和ループを通さず既存式をそのまま実行する分岐にする —
     `x_i=1.0` の乗算 1 回でも浮動小数点上は同値だが、演算列を変えない
     方が安全)。
4. **(Z1,Z2) ペア別散乱テーブルキャッシュ** (確定事項):
   ```cpp
   // Cache of per-(dopant Z1, target Z2) scatter tables. Built single-threaded
   // in prebuild_scatter_tables() BEFORE the worker loop; read-only afterwards.
   std::map<std::pair<int,int>, ScatterTable>& scatter_cache();
   const ScatterTable& get_scatter_table(int z1, int z2);  // insert if absent
   ```
   現行の ZBL 万能遮蔽では全エントリは数値的に同一だが、キャッシュ構造を
   入れておく (将来のペア別ポテンシャル差し替え点。overview の確定構造)。
   `apply_mc_implant` 内 L465 の `scatter_table();` を
   `for (auto& mc_ : w.mats) for (auto& c : mc_.comps) c.table =
   &get_scatter_table(z1, c.z);` に置換 — **worker 起動前にメイン
   スレッドで一括構築** (map への挿入は並列域外のみ。並列域内は
   ポインタ経由の読み取り専用でロック不要)。既存の
   `mc::sin2_half_theta` (L166) は万能テーブルを引く現行実装のまま残す
   (test_mc.cpp の検証 API)。
5. `walk_ion` の衝突相手抽選 (L370–375 を置換):
   ```cpp
   // Collision partner: proportional to number fraction x_i.
   // CRITICAL (determinism): the partner draw consumes one rng.u() ONLY
   // when comps.size() > 1. Single-element materials must take the
   // comps.size()==1 branch with ZERO extra draws so the per-ion RNG
   // stream is bit-identical to the current implementation.
   const CompConstants* cc = &mc_.comps[0];
   if (mc_.comps.size() > 1) {
     const double u = rng.u();               // partner draw
     for (const auto& c : mc_.comps)
       if (u < c.x_cum) { cc = &c; break; }  // x_cum precomputed cumsum
   }
   const double p_sq = b_min_sq + rng.u()*(mc_.pmax*mc_.pmax - b_min_sq);
   const double p    = std::sqrt(p_sq);
   const double s2   = cc->table->sample(e * cc->eps_per_ev, p * cc->inv_a);
   const double t_recoil = cc->tmax_fac * e * s2;
   ```
   **RNG ストリーム内の位置を仕様として固定する**: 相手抽選は
   「衝突径数 draw の直前」。1 イオンの draw 列は
   `pos.x, pos.y, {[partner], impact, azimuth}×ステップ` であり、
   `[partner]` は多成分材料内のステップのみ存在する。散乱角計算 L389–392 の
   `mass_ratio` も `cc->mass_ratio` に差し替え。
6. ダメージ (L379–387) は現行どおり `crystal_si` 材料のみ・閾値 `kEdSi`
   のまま (化合物材料は crystal_si=false なのでダメージ蓄積なし —
   現行挙動を変えない)。チャネリング (L347–365) も crystal_si のみで
   無変更。
7. 新しい調整可能定数は導入しない。もし較正の必要が出た成分比・変位
   閾値等は P1-10 の `ParamDB` (`include/cprocess/param_db.hpp`,
   `ParamDB::instance().get("mc.<key>", fallback)`) キーとして読むこと。
   本タスクの範囲では ParamDB キー追加なし。

## テスト仕様

新規 `tests/test_compound_bca.cpp` (main + `test_util.hpp` の assert
パターン)。`CMakeLists.txt` L76 の `foreach(t ...)` に `compound_bca` を追加。
共通設定: `channeling=false` (決定性モード)、`threads=1` と `4` の両方で
seed 固定。メッシュは `make_box_mesh(0,0.2e-4, 0,0.2e-4, 0,0.5e-4, 4,4,50)`
程度、`silicon_mask` は全 1。

1. **Si ビット一致 (退化ケース)**: 同一パラメータ (B, 30 keV, dose 1e14,
   ions 50000, seed 42) で (a) `material_table = {target_silicon()}` +
   cell_material 全 0、(b) comp を明示した
   `TargetMaterial t = target_silicon(); t.comp = {{14, 28.086, 1.0}};`
   の 2 通りを実行し、`conc` 配列の全要素と `McImplantStats` の全整数
   カウンタ・`rp`/`drp` が **ビット一致** (`==` 厳密比較) であること。
   さらに (c) cell_material=nullptr のレガシーパス (L433–434) とも
   ビット一致であること。threads=1 と threads=4 で同一。
2. **SiO₂ 中 B 30 keV の Rp**: 全セル oxide の cell_material で
   (a) 有効単一元素版 (`comp` を空にした target_oxide()) と
   (b) 化合物版 (成分 {Si 1/3, O 2/3}) を ions=200000, seed=7 で実行。
   `rp_a`, `rp_b` について
   `0.05 <= |rp_b - rp_a| / rp_a <= 0.20` かつ
   `|rp_b - 100e-7| < |rp_a - 100e-7|` (SRIM 参照値 ~100 nm へ近づく方向)
   を assert。
3. **軽元素後方散乱の増加 (定性)**: B 5 keV, ions=200000 を
   (a) Si 単体、(b) SiO₂ 化合物で実行し、
   `stats_b.backscattered > stats_a.backscattered` を assert
   (O 反跳による大角散乱の寄与)。
4. **実行時間 < 30% 増**: 物理的に同一な 2 ケース
   (a) Si 成分 1 個、(b) Si を 2 成分 `{{14,28.086,0.5},{14,28.086,0.5}}`
   に分割 (抽選パスのオーバーヘッドのみが差) を ions=500000, threads=1 で
   `std::chrono::steady_clock` 計測し、`t_b < 1.3 * t_a` を assert。
   マシン負荷ノイズ対策として 3 回計測の最小値同士を比較する。
5. **既存回帰**: `test_mc` / `test_resist` / `test_photo` /
   `test_integration` を含む既存全テストが無変更で PASS すること
   (単一元素パスのビット一致が保証するので既存テストの期待値更新は不可)。

Python 側: `proc::` シグネチャ変更なしのため新規バインディング不要。
`python/test_comprehensive.py` には photo()+MC 実装の既存テストが SiO₂ を
使わないことを確認の上、追加不要 (mat_table に oxide を積む経路が
できたら別タスク)。

## 完了条件 (DoD)

- [ ] `TargetComponent` / `TargetMaterial::comp` が mc_implant.hpp に追加、
      既存の集約初期化がコンパイル互換
- [ ] 衝突相手抽選が x_i 比例、成分 1 個では RNG 追加消費ゼロ
- [ ] (Z1,Z2) ペア別 ScatterTable キャッシュが worker 起動前に一括構築
- [ ] `target_oxide()` が Si/O 2 成分、`target_nitride()` 新設
- [ ] 電子的阻止能が Bragg 加重和、自由行程が総数密度
- [ ] `tests/test_compound_bca.cpp` の 4 テスト PASS、CMakeLists 登録済み
- [ ] 既存全テスト PASS (`ctest --test-dir build`)、特に test_mc は
      期待値無変更で PASS (ビット一致の証明)
- [ ] OpenMP 無効ビルドでもコンパイル・PASS
- [ ] コミットメッセージに `P3-a` を含める

## やらないこと

- 反跳カスケードの追跡 (recoil を二次イオンとして輸送しない —
  ダメージは現行 Kinchin-Pease カウントのまま)
- 元素別変位閾値 E_d (kEdSi=15eV のまま; 化合物は crystal_si=false で
  ダメージ対象外)
- 化合物中のチャネリング (crystal Si のみ、現行どおり)
- ZBL 以外のペア別ポテンシャル (キャッシュ構造だけ用意、中身は万能 ZBL)
- `proc::` / pybind / Simulation のシグネチャ変更 (材料表は内部データ)
- スパッタリング・表面ミキシング
