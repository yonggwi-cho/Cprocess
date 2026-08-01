# C-1: 注入モーメント表の拡充 + dual-Pearson

## 目的

解析注入 (`proc::implant_gauss`) の Pearson-IV モーメント表は 8 種×5-8 点
(10-200 keV) しかカバーせず、現代的な低エネルギー(源/ドレイン extension、
~1 keV)から高エネルギー(ウェル/リトログレード、~MeV)実装の大半が範囲外
クランプで潰れる。加えて単峰 Pearson-IV は MC(BCA)が示す
"チャネリングテール"(結晶軸に沿った低角散乱イオンが Rp の数倍先まで
残す長い裾)を全く再現できず、`tests/test_sprocess_parity.cpp` の
`[C-1] analytic channeling tail @2Rp` が解析値を MC の 1e-7〜1e-8 倍
過小評価していた(v2 分析 物理 Top1 ギャップ)。本タスクはモーメント表を
1 keV〜3 MeV・15 点以上/種に拡充し、`profile="dual"`(主峰 Pearson-IV +
指数チャネリングテール)を追加してこのギャップを埋める。

## 現状コード

- `include/cprocess/materials.hpp` — `Dopant::range`
  (`std::vector<std::array<double,5>>` = {E keV, Rp cm, dRp cm, γ, β});
  `implant_moments(d, E, rp, drp, gamma, beta)` が log-E 線形補間、範囲外
  クランプ
- `src/materials.cpp` — `kDopants` テーブル。B/P/As/C/F: 8 点
  (10-200 keV)、Sb/In/Ge: 5-8 点。B/P/As/Sb/In/C/F/Ge の 8 種
  (C/F/Ge も range テーブルを持つ — 「C/F は range テーブル不要」という
  前提は誤りで、実際には既に全種が保持している)
- `include/cprocess/implant.hpp` — `ImplantParams{ profile: gauss|pearson4 }`
- `src/implant.cpp` — `build_pearson4()`(Pearson-IV ODE 積分)、
  `apply_implant()`(profile==pearson4 ならテーブル評価、それ以外は
  Gaussian)
- `src/process.cpp` — `proc::implant_gauss(..., profile="gauss"|"pearson")`;
  `profile=="pearson"` で `implant_moments()` を引いて
  `ImplantParams::Profile::pearson4` を組み立てる
- `tests/test_sprocess_parity.cpp` — 失敗中の `[C-1]`
  チェック(B 40 keV, `implant_mc(channeling=true)` を基準に 2·Rp での
  濃度比を検証)

## 実装手順

### (1) モーメント表拡充(`src/materials.cpp`)

8 種全てについて、既存の 10-200 keV アンカー行(Tier-A ベンチマークの
恒久回帰アンカーであり**ビット不変で温存**)はそのまま残し、
1/2/5 keV(低エネルギー側)と 300/400/600/800/1000/1500/2000/2500/3000 keV
(高エネルギー側)を追加して 17〜20 点/種にした。

**データ導出方法**(コード内コメント `src/materials.cpp` 冒頭にも記載):

- 電子阻止能側の漸近: `src/mc_implant.cpp` の `kls_coeff` が実装している
  Lindhard-Scharff 電子阻止能 dE/dx ∝ √E と同じ物理から、
  高エネルギー極限で射影飛程 R(E) ∝ √E(対数-対数勾配 n=0.5)。
- **高エネルギー側 (>200 keV, 〜3 MeV)**: 表の最後の 2 点(150/200 keV)
  から測った局所対数-対数勾配 n を、ln(E) に関して 3 MeV で n=0.5 に
  線形に緩和させ、各新点を段階的べき乗則
  `Rp(E_i) = Rp(E_{i-1})·(E_i/E_{i-1})^n(E_i)` で計算(dRp も同様)。
- **低エネルギー側 (<10 keV, 1 keV まで)**: 表の最低 2 点間の局所勾配
  (核阻止能が飽和し概ね線形に近い領域)を 1 keV まで一定に保って外挿
  (核阻止能の Lindhard 普遍関数は reduced energy ~0.3-1 に幅広い極大を
  持ち、この付近では R がほぼ E に比例し続ける)。
- **γ/β(歪度・尖度)**: 上記のような厳密な閉形式 LSS モーメント式は
  本コードベースの他所にも存在しないため、**近似**として扱う。
  高エネルギー側は 3 MeV に向けて Gaussian 極限 (γ→0, β→3) に緩和
  (高エネルギーほど核散乱由来の非対称性が薄れるという定性的傾向は
  既存の 200→10 keV 行でも見られる)。低エネルギー側は既存の
  10/20 keV 2 点間の ln(E) 線形傾向を延長し、γ∈[-3,0]・β∈[3,12] に
  クランプ。
- 既存の 10-200 keV コア値はどの新点によっても**変更していない**
  (`implant_moments` の線形補間はテーブル値そのものを通すため、
  この区間のクエリはビット不変)。

### (2) dual-Pearson プロファイル

- `include/cprocess/implant.hpp`: `ImplantParams::Profile` に `dual` を
  追加。新フィールド `dp_frac`(テールへのドーズ分配率)・`dp_l`
  (テールの指数減衰長 [cm])。
- `src/implant.cpp` の `apply_implant()`: `profile==dual` のとき、
  主峰 Pearson-IV(ドーズの `(1-dp_frac)`)に加え、`d >= rp` の領域で
  ```
  C_tail(d) = dp_frac * dose / dp_l * exp(-(d - rp) / dp_l)
  ```
  を加算する一枚指数のチャネリングテール(片側、`d<rp` では 0)。
  `dp_frac<=0` または `dp_l<=0` の場合はテールが消え、通常の Pearson-IV
  (もしくは Type-IV 無効時の Gaussian フォールバック)にビット一致で
  縮退する。
- `src/process.cpp` の `implant_gauss()`: `profile=="dual"` を受理し、
  `implant_moments()` で主峰モーメントを引いた上で、既定値を
  ParamDB キー `<Sym>.dp.frac`(既定 0.06)・`<Sym>.dp.decay_mult`
  (既定 2.0、`dp_l = decay_mult * drp`)から読む。Python 側は
  `profile=` 文字列パラメータをそのまま `_c.proc_implant_gauss` へ
  渡す既存の仕組みで `"dual"` も透過するため、シグネチャ変更不要
  (CLAUDE.md のユニット変換ルールは元々 `implant()` が満たしている)。

### (3) 較正

`proc::implant_mc(channeling=true)` を較正オラクルとして使用。B の
5/20/40/80 keV について解析(dual, `dp_frac=0.06`, `dp_l=2·dRp`)と MC
(200k〜400k イオン)の深さ方向濃度を比較(スクラッチスクリプト、
コミットせず、この文書に手順を記載——再現可能):

| E [keV] | Rp(解析) [nm] | 深さ倍率 | 解析濃度 [cm⁻³] | MC濃度 [cm⁻³] | 比 |
|---|---|---|---|---|---|
| 40 | 148.7 | 2.00×Rp | 1.10e17 | 8.91e17〜9.99e17(seed依存) | 0.11〜0.12 |
| 20 | 73.75 | 1.75×Rp |  | 8.72e18 | (定性的に同様の減衰形) |
| 80 | 271.3 | 1.50×Rp |  | 4.24e18 | (同上) |

`dp_frac=0.06`・`dp_l=2·dRp` は B 40 keV での parity チェック
(2·Rp での比 ≥ 0.10)を安全マージン込みで満たす(実測比 0.11-0.12、
要求 0.10)。他エネルギーでの厳密な L2 一致(<15%)は将来の較正課題として
残る(dp_frac/dp_l は ParamDB で種別に上書き可能なので、他種・他エネルギー
の較正はロジック変更なしで追加できる)。

## テスト仕様

- `tests/test_dual_pearson.cpp`(CTest 名 `dual_pearson`、通常スイート・
  常時 PASS):
  1. **[C-1] parity チェック本体**(`test_sprocess_parity.cpp` から移設):
     B 40 keV, `profile="dual"` vs `implant_mc(channeling=true)`、
     2·Rp での濃度比 ≥ 0.1。
  2. `profile="dual"` のドーズ保存(±1%)。
  3. `dp_frac=dp_l=0`(既定 `ImplantParams`)なら `pearson4` と
     ビット一致(フォールバック不変量)。
  4. ParamDB `<Sym>.dp.frac`/`<Sym>.dp.decay_mult` の上書き可能性:
     `frac=0` でテール消失、`frac=0.20` でテール増大を確認。
  5. 拡充モーメント表の定性チェック: 1 keV〜3 MeV で Rp が単調増加、
     3 MeV 近傍の対数-対数勾配が LSS 高エネルギー漸近値 0.5 近傍
     (0.3〜0.7)、3 MeV での γ/β が Gaussian 極限近傍、1 keV の Rp が
     10 keV アンカーより小さいこと。
- 既存 `tests/test_pearson.cpp`・`tests/test_benchmarks.cpp` の
  Tier A 12 項目(B/P/As の Rp/dRp × 30/50/100 keV)は無変更で再実行し
  全 PASS を確認(10-200 keV アンカー行がビット不変のため回帰なし)。

## 完了条件 (DoD)

- [x] `implant_moments()` が全 8 種で 1 keV〜3 MeV・17 点以上を返す
- [x] `profile="dual"` が実装され、B 40 keV で parity チェック
      (2·Rp で解析 ≥ MC/10)が PASS
- [x] `tests/test_dual_pearson.cpp` が `CMakeLists.txt` に登録され PASS
- [x] `test_pearson`/`test_benchmarks`(Tier A 全 12 項目+他)が無回帰で PASS
- [x] `tests/test_sprocess_parity.cpp` から `[C-1]` チェックを除去し、
      `docs/tasks/README.md` のパリティ表・要約行を更新
- [x] python バインディング/`Simulation` 側はシグネチャ変更不要
      (既存の `profile=` 文字列引数がそのまま "dual" を通す)ため
      追加実装なし——ただし `Simulation.implant()` の docstring に
      `profile="dual"` の説明を追記

## エンジニアはこの工程の結果をどう確認するか

- **ログ**: `implant_gauss()` のログ行が
  `profile=dual(pearson4+tail)`(または Type-IV 無効時
  `dual(gauss+tail)`)を出力し、主峰形状とテールの有無が一目でわかる。
- **抽出**: `Simulation.profile1d()`(既存)で深さ方向濃度プロファイルを
  取得し、`profile="pearson"` と `profile="dual"` を重ねてプロットすれば
  Rp 近傍は同一、Rp の 1.2〜2 倍超で dual 側だけ緩やかなテールが立つ
  ことが目視確認できる。
- **可視化**: 既存の `Simulation.save()`/VTU 出力に濃度場がそのまま
  乗るため、ParaView 等での断面表示でも確認可能(dual はテール分だけ
  裾野が沈み込むオレンジ~赤の帯として見える)。
- **パラメータ確認**: `ParamDB.all()`(既存)で `<Sym>.dp.frac`/
  `<Sym>.dp.decay_mult` の現在値(既定 or 上書き済み)を確認できる。

## やらないこと

- SRIM の実行(サンドボックスで不可能) — 既存コードベースの LSS/
  Lindhard-Scharff 物理(`kls_coeff` 前例)からの解析的外挿に限定
- γ/β の厳密な LSS/BCA 由来モーメント式の実装 — Gaussian 極限への
  緩和という定性的近似に留める(コード内で明記)
- dual テールの tilt/rotation 依存(横方向チャネリングは対象外、縦方向
  深さプロファイルのみ)
- B 以外の種(P/As/Sb/In/Ge/C/F)での MC 較正 — ロジックと ParamDB
  キーは種非依存で用意したが、既定値のフィットは B のみ実施
  (他種は `<Sym>.dp.*` で個別上書きして追って較正)
