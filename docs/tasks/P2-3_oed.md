# P2-3: OED/ORD — 酸化増速拡散 (酸化界面からの I 注入)

## 目的

酸化中の Si/SiO₂ 界面は成長速度に比例して格子間原子 (I) を Si 側へ注入し、
B/P の拡散を増速する (OED)。Sb など空孔支配種は逆に遅延する (ORD —
I 過剰が I-V 再結合で V を減らすため、P2-1 の連立モデルで自動的に出る)。
本タスクで `proc::oxidize` に I 注入源を追加し、点欠陥緩和と結線する。

**P1-6 (oxidize 統合) と P2-1 (フル点欠陥) の完了が前提。**

## 現状コード

- `src/process.cpp` — `proc::oxidize(st, time_s, temp_k, wet, log)` (P1-6):
  Deal-Grove で膜厚増分 `dx_ox` を計算し、`extend_mesh_exact` + 再タグで
  形状実現。**時間分割はしていない** (1 回の幾何更新)。
- `src/diffusion.cpp` — `run_ted` (P2-1 後): `st.fields["I"]` (過剰濃度) を
  初期値として C_I を解く。`proc::diffuse_ted` が "I"/"V"/"C311" を結線済み。
- `src/oxidation.cpp` — `deal_grove_step(x0_um, dt_min, T_c, wet)`。
- `src/materials.cpp` — Si 原子密度は mc_implant.cpp の `kNt = 4.99e22`
  にあるのみ。本タスクで `constexpr double kNSi = 5.0e22;` を
  materials.hpp に公開する。
- P1-10 ParamDB: `params.get("oed.theta", 0.01)`。

## 実装手順

1. **離散 I 注入源の定義**。1 サブステップで酸化膜が `d(dx)` 成長したとき、
   界面から注入される I の面密度は `theta · d(dx) · N_Si` [cm⁻²]
   (theta = 注入率、既定 0.01、ParamDB キー `"oed.theta"`)。これを
   「Si/SiO₂ 界面直下 1 セル層」の Si セルに体積ソースとして与える:
   ```cpp
   // 界面高さ z_if (P1-6 が計算済み) の直下、z_if - h_cell <= centroid.z < z_if
   // にある silicon セル i に対し:
   //   fields["I"][i] += theta * d_dx * kNSi / h_cell;   // cm^-3
   // h_cell = infer_box_dims の z セル高 (base_lz / nz)
   ```
   セル層の判定はセル centroid で行う。マスク開口 (将来 P2-4) に備え
   x/y は制限しない (ブランケット)。この式をコードコメントに明記する。
2. **oxidize のサブステップ化 + 点欠陥緩和の内蔵**。
   `proc::oxidize` を次のループに変更 (N = 10 固定サブステップ):
   ```
   for k in 1..N:
     dt_k = time_s / N
     (a) Deal-Grove: x_{k} = deal_grove_step(x_{k-1}, dt_k, ...) ; d_dx = x_k − x_{k-1}
     (b) 幾何更新: d_dx 分の extend/再タグ (P1-6 の (d)〜(h) を増分適用。
         rise が z セル高の 1/4 未満のサブステップは幾何更新を繰越し累積
         (pending_rise) し、閾値を超えたステップでまとめて実現する —
         メッシュを 10 回作り直さないための実務上の措置)
     (c) I 注入: 手順 1 の式で fields["I"] に加算 (幾何繰越しに関係なく毎回)
     (d) 点欠陥+ドーパント緩和: run_ted 相当を dt_k だけ実行。実装は
         proc::diffuse_ted(st, opts_k) を呼ぶ (opts_k: temp=temp_k,
         time=dt_k, verbosity=0)。ドーパント場が無ければ点欠陥のみ進む。
   ```
   これにより **単一の oxidize() 呼び出しだけで OED が現れる**。
   また oxidize が "I" 場を書くので、後続の `diffuse_ted` も自動で拾う。
3. **theta=0 の分岐**: `oed.theta == 0` なら (c)(d) をスキップし
   P1-6 の従来動作 (幾何のみ) と一致させる。ログに
   `[oxidize] OED theta=... I_injected=...` を追加。
4. **proc:: シグネチャは不変** (oxidize の引数追加なし; theta は ParamDB)。
   よって pybind/Simulation の変更は不要。Python テストのみ追加。

## テスト仕様

`tests/test_oed.cpp` (新規、CMakeLists foreach に `oed` 追加)。
メッシュ: `mesh_box(st, 0, 0.2e-4, 0, 0.2e-4, 0, 0.5e-4, 4, 4, 50)`、
B 1e18 cm⁻³ を深さ 0.1 µm に Gaussian 埋込み (implant_gauss rp=0.1e-4)。

1. **OED 増速**: state A: `oxidize(30 min, 1000 °C, wet)` /
   state B (inert): `set_param("oed.theta", 0)` で同じ oxidize
   (幾何+熱履歴同一)。両者の B の深さ方向 spread (test_ted.cpp の
   `profile_spread` を流用) を比較し **spread_A > 1.2 × spread_B**。
2. **theta=0 で inert と一致**: state B と、oxidize せず
   `diffuse_ted(30 min, 1000 °C)` のみの state C の spread が相対差 < 2%
   (幾何変化分のセル転写誤差を除くため、spread は Si セルのみで評価)。
3. **dose 保存**: state A の B 総量 (Σ C·V、oxide 転化セル含む) が
   oxidize 前後で相対差 < 5% (P1-6 の転写許容と同じ)。
4. **I 場が正**: state A の oxidize 直後 `max(fields["I"]) > 0` かつ
   界面直下セルで最大 (argmax セルの centroid.z が z_if − 2h 以上)。

Python (`python/test_comprehensive.py`):
`sim.init("B", 1e18).oxidize(30, 1000, wet=True).diffuse(5, 1000, ted=True)`
と、`sim.set_param("oed.theta", 0.0)` した同一フローの B spread を numpy で
比較し前者が大きいこと。既存 oxidize テストが PASS すること。

## 完了条件 (DoD)

- [ ] oxidize が N=10 サブステップで 成長→注入→点欠陥緩和 を回す
- [ ] 注入式が仕様どおり (`theta·d_dx·N_Si/h_cell`、ParamDB `oed.theta`)
- [ ] theta=0 で P1-6 従来動作と一致 (テスト 2)
- [ ] C++ テスト 4 本 + Python テスト PASS、CMakeLists 登録、全テスト PASS
- [ ] コミットメッセージに `P2-3` を含める

## やらないこと

- V 注入 (ORD は P2-1 の I-V 再結合経由で間接的に出るもののみ)
- マスク開口酸化での横方向 OED 分布 (P2-4 の 2D 酸化後に自然拡張)
- theta の温度/雰囲気依存 (定数 1 個。較正は ParamDB で)
- 偏析・dose loss の変更 (P1-4 のまま)
- oxidize の引数追加 (theta はパラメータ機構経由のみ)
