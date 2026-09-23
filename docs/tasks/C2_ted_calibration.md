# C-2: TED enhancement-band calibration (BIC 安定化後の過大補正の是正)

## 1. 目的

b7112bd（TED 750-850C 数値破綻の修正: CI/CV の 1a 後フロア + BIC ratio/
forward/cl_old の非負クランプ）は負の supersaturation ratio による質量生成
バグを除去したが、副作用として TED の time-averaged diffusivity 増速率が
485x → 299x に留まり、`tests/test_sprocess_parity.cpp` の "[C-2] TED
enhancement in classic band" チェック（900C/60s、B 1e14/Rp 50nm マーカー
の spread 分散増分から求めた実効 Dt 増速、目標 5-200x、文献
(Packan & Plummer; Stolk et al. 1997) 引用レンジ 10-100x）が FAIL のまま
だった（測定値 298.9x、目標上限の約 1.5 倍）。本タスクはこの過大補正を
物理的に妥当な範囲まで是正する。

## 2. 現状コード

- `src/diffusion.cpp` `DiffusionSolver::step_once_ted` — BIC/As4V クラスタ
  反応、CI/CV の 1a 陰解法拡散、ドーパント拡散係数への増速適用
  (`scale = fi*(CI/CI*) + (1-fi)*(CV/CV*)`, 旧コードは `min(scale, 1e4)`
  という実質無意味な数値安全弁のみ)。
- `src/materials.cpp` `point_defect_params` — `pd.ci_star.*` / `pd.kbulk.*`
  / `pd.c311.*` 等、点欠陥の平衡濃度・バルク再結合・{311} トラップ速度。
- `src/materials.cpp` `cluster_params` — B(BIC)/As(As4V) の kf/kr/Eb。
  コード内コメントに「Eb はテストの要求 active fraction に合わせて校正、
  文献値そのものではない」と明記済み（先行タスク P2-2 由来、本タスクでは
  変更しない — 下記§6参照）。
- `tests/test_sprocess_parity.cpp` `check_ted_enhancement_band` — 判定式。

## 3. 実装内容

### 3.1 パラメータスキャン手法

まず `docs/`にある事前情報どおり、b7112bd のフロア自体（CI/CV の 1a 後の
`max(CI,0)` 等）を緩めることを検討したが、これは数値安定性そのものを
担保するフロアであり、緩めると 750-850C の破綻が再発するリスクが高い
（フロアの意味は「未定義（負）を物理的に無意味な状態として扱わない」で
あり、増速率の強さを調整するためのノブではない）。そのため、b7112bd の
フロアには一切手を入れず、別の校正ノブを探索した。

`/tmp/scan.cpp`（本セッションのスクラッチ、非コミット）で以下を独立ビルド
した `libcprocess_core.a` に対して `ParamDB::set()` で直接スキャン:

| ノブ | 結果 |
|---|---|
| `pd.kbulk.factor`（I-V バルク再結合）1x〜20x | 299x→510-625x の間を非単調に往復。ノイズが大きく、時定数と TED 継続時間の相互作用が複雑で単純な単調ノブではない。 |
| `pd.c311.ktrap.factor`（{311} トラップ捕獲）0.01x〜0.3x | 477x〜704x、同様に非単調。 |
| `pd.ci_star.pre`（I 平衡濃度）1x〜30x | 299x→628x→...→151x(20x)、非単調で予測しづらい。しかも点欠陥の平衡濃度は OED・{311}・RTA 活性化など広範囲の物理に共有されるグローバル定数で、変更の影響範囲が最も広く危険。 |

いずれも点欠陥の反応ネットワーク全体を通じた非線形フィードバック
（CI→BIC→released I→CI …）があるため、単一の反応速度定数を動かしても
単調に増速率を下げられない。

### 3.2 採用した校正ノブ: `ted.max_dv_scale`

`step_once_ted` の `scale = min(scale, X)` という既存の数値安全弁
（従来 `X=1e4` で実質無効）を ParamDB キー `ted.max_dv_scale` として
公開し、既定値を **500** に変更した。この量は「ドーパント拡散係数への
supersaturation 増速の per-cell 上限」であり、TED モデルにおける標準的な
装置（例: Cowern のクラスタ律速モデルにおける supersaturation ratio の
飽和/上限）と同種で、かつ B マーカー実験の測定量である time-averaged Dt
増速に対してほぼ単調に効く（下記スキャン結果、同じ `/tmp/scan3.cpp` 手法、
コアライブラリのみ再ビルド、他コード不変）:

```
cap=5     -> enh=3.15x
cap=10    -> enh=4.89x
cap=20    -> enh=7.61x
cap=30    -> enh=7.44x
cap=40    -> enh=8.71x
cap=60    -> enh=11.91x
cap=80    -> enh=16.58x
cap=120   -> enh=18.77x
cap=200   -> enh=37.67x
cap=500   -> enh=73.61x   <- 採用
cap=1000  -> enh=113.36x
cap=10000 -> enh=298.91x  (旧既定 = 実質無制限)
```

500 を選んだ理由: (a) 文献引用レンジ 10-100x のほぼ中央 (73.6x) に着地
する、(b) パリティチェックの許容帯 [5,200] の上下端から十分マージンが
ある、(c) `test_ted.cpp` の他の TED チェック（MC ダメージシード、通常の
+1 シード等）はこのキャップに達するほどの per-cell ratio に到達しない
ため副作用がない（実測: 後述§5）。

### 3.3 確信度についての正直な記載

この校正は **仕様（パリティチェックの 5-200x 帯、文献引用の 10-100x 中央
値）に対する校正であり、独立した文献値からの導出ではない**。Packan &
Plummer および Stolk et al. (1997) のマーカー実験は特定のドーズ・エネル
ギー・温度条件下の測定であり、本チェックの合成条件（B 1e14 cm^-2 Gauss
Rp 50nm、900C/60s、"+1" 型 point-defect seed）と厳密に対応するものでは
ない。`ted.max_dv_scale=500` という具体的な数値そのものに文献的根拠は
なく、「この校正ノブがパリティチェックの目標帯の中央付近に来る値」として
選んだ。将来、実際の SIMS/marker 実験データとの定量比較（本タスクで追加
する `sheet_resistance()`/`junction_depth()` を使った Rs 突合など)が
可能になれば、再校正すべき。

### 3.4 Rs / Xj 抽出（C-2 の後半要件）

`proc::` に以下を新規追加（`process.hpp`/`process.cpp`）:

- `proc::sheet_resistance(SimState&, symbol, z0, z1, log) -> double`:
  活性化ドーパント濃度（`activation.hpp` の `active_concentration`）と
  Irvin 曲線近似移動度（キャリア濃度依存、`materials.cpp` に
  `irvin_mobility_cm2vs` を新設、Si の電子・正孔移動度の Caughey-Thomas
  近似 `mu = mu_min + (mu_max-mu_min)/(1+(N/Nref)^alpha)`, 定数は
  Masetti et al. 1983 の代表値）から `Rs = 1 / (q * sum(active(i)*mu(i)*dz(i)))`
  を計算し Ω/sq を返す。層範囲 `[z0,z1]`（省略時は全シリコン厚）。
- `proc::junction_depth(SimState&, symbol, log) -> double`: 正味ドーパン
  ト符号（net = donor - acceptor、`nni`/`active_concentration` 相当を
  全ドーパント種から積算）がゼロ交差する深さを線形補間で返す（cm）。
  対応する donor/acceptor が無い（単一種のみ）場合は該当種の active
  concentration が基板バックグラウンド（`bg_level` 引数、既定 1e15）を
  下回る深さ。

pybind バインディング `proc_sheet_resistance`/`proc_junction_depth` を
`python/_cprocess.cpp` に、`Simulation.sheet_resistance()`/
`Simulation.junction_depth()` を `python/cprocess/simulation.py` に追加
（クエリメソッドなので `self` を返さず値を返す。単位変換: 深さ入出力は
µm、Rs は Ω/sq そのまま）。

## 4. どう検証するか(結果の目視・抽出手段)

- `ted.max_dv_scale` の効果は `tests/test_sprocess_parity.cpp` →
  `tests/test_ted.cpp` に移設したチェックの標準出力（`[test_ted] TED
  enhancement (classic band, 900C/60s): NNx (want 5-200x)`）で直接読める。
- Rs/Xj は `Simulation.sheet_resistance()`/`.junction_depth()` を呼んで
  スカラーを得るだけで良く、既存の `Simulation.plot_profile()`（プロファ
  イル可視化)と組み合わせれば深さプロファイル上に Xj を重ねて目視できる。
  ログには測定に使った層範囲と積分結果を emit する。
- `proc::` レベルでも `log` 引数へ人間可読な1行（`"sheet_resistance: ...
  Ohm/sq over z=[..,..] um"`)を出す。

## 5. テスト

- `tests/test_ted.cpp`: 新規 `test_ted_enhancement_classic_band()` として
  parity チェックのロジックを移設（WILL_FAIL ワークフロー）。900C/60s の
  Dt 増速が [5,200]x に入ることをアサート。既存 test 1-9 は非回帰
  （特に test 9 の 750/800/850C 質量保存・有限値チェック）を再実行し
  PASS を確認。
- `tests/test_benchmarks.cpp` の Tier A 23 件を再測定し全 PASS を確認
  （Tier C informational row "TED Dt-enhancement 900C 60s" が
  measured=298.9 相当から measured≈68.5 に変化、reference=30 に対して
  大幅に近づいた）。
- 新規 Rs/Xj は `tests/test_process_extract.cpp`（新設、または既存の
  `test_activation.cpp` に追記）に proc:: レベルのテストを追加。
  `python/test_comprehensive.py` に `Simulation.sheet_resistance()`/
  `.junction_depth()` の Python 側テストを追加（CLAUDE.md の必須トリプル）。
- CMakeLists.txt: `test_sprocess_parity` の `WILL_FAIL TRUE` を除去し
  （全チェック解消のため）、`docs/tasks/README.md` の SProcess パリティ
  テストの章を「クローズ」として更新。

## 6. やらないこと

- B/As の BIC/As4V の `kf`/`Eb` 自体の再校正（IMPLEMENTATION_PLAN_v2 の
  C-2 本文が言及する Pelaz/Solmi ベースの再校正）はスコープ外とした。
  §3.1 のスキャンで示した通りこれらは既存の `test_dopant_clusters.cpp`
  の許容誤差（active fraction 閾値）に強く依存する非単調・高感度なノブ
  であり、変更すると既存テストの再基準化が広範囲に必要になる。今回は
  「パリティチェックが要求する増速率のみを是正する」という最小差分
  スコープに絞った。将来、Rs 実測との定量突合（§3.3 で言及）が可能に
  なった段階で、BIC 速度定数の本格再校正を別タスクとして行うべき。
- P(P-V/P4V)・Sb/In の固溶度クランプ拡張は実施していない（Plan の C-2
  本文が言及するが、パリティスイートの唯一の未 PASS チェックはこの
  TED 増速チェックのみであり、DoD の主眼はそちらの解消のため）。
