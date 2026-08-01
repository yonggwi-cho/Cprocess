# C-3: 酸化の薄膜補正・雰囲気依存 + Deal-Grove 係数の ParamDB 化

## 1. 目的

`proc::oxidize()` の Deal-Grove 係数はこれまで `src/oxidation.cpp` にハード
コードされており、(a) 実験較正やパラメータスタディでの調整ができない、
(b) 10 nm 未満の薄膜域で実測より系統的に薄く成長する(Massoud 補正の欠如)、
(c) 圧力・HCl 添加・結晶方位といった SProcess が標準サポートする雰囲気依存
性が一切モデル化されていない、という 3 つのギャップがあった
(`docs/IMPLEMENTATION_PLAN_v2.md` C-3、`tests/test_sprocess_parity.cpp` の
"Massoud thin-oxide enhancement" WILL_FAIL チェック)。

本タスクは (1) Deal-Grove の B/A アレニウス定数を ParamDB 化(既定値は現行
定数と厳密一致、ビット不変のリファクタ)、(2) Massoud (1985) 型の薄膜増速項
を ParamDB キー付きで追加(既定 OFF、opt-in でビット不変性を壊さない)、
(3) 圧力スケーリング・HCl 増速・方位係数を `proc::oxidize()` の新規オプショ
ン引数として追加する。

## 2. 現状コード

- `src/oxidation.cpp` — `deal_grove_dry`/`deal_grove_wet`: B, B/A の
  アレニウス前指数・活性化エネルギーがハードコード。`deal_grove_step`:
  閉形式積分(1ステップ厳密解)。
- `include/cprocess/oxidation.hpp` — 上記のシグネチャ。
- `src/process.cpp` — `proc::oxidize()`(1D ブランケット酸化、P1-6/P2-3):
  `theta<=0` の場合は `deal_grove_step` を1回、`theta>0` の場合(既定)は
  N=10 サブステップに分けて `deal_grove_step` を繰り返し呼ぶ(P2-3 OED)。
  いずれも Deal-Grove 係数は方位・圧力・HCl 非依存の <100>/1atm/HCl なし
  固定だった。
- `include/cprocess/process.hpp` — `oxidize()` 宣言(P1-6/P2-3 のコメント
  参照)。
- `python/_cprocess.cpp` の `proc_oxidize` バインディングと
  `python/cprocess/simulation.py` の `Simulation.oxidize()`。

## 3. 実装内容

### 3.1 Deal-Grove 係数の ParamDB 化

`deal_grove_dry`/`deal_grove_wet` の前指数・活性化エネルギーを ParamDB
キーから読むよう変更(`ParamDB::instance().get(key, fallback)`、fallback
は現行ハードコード値そのもの):

| キー | 既定値 | 意味 |
|---|---|---|
| `ox.dry.b0` | 7.72e2 | 乾燥酸化 B 前指数 (um²/hr) |
| `ox.dry.be` | 1.23 | 乾燥酸化 B 活性化エネルギー (eV) |
| `ox.dry.a0` | 6.23e6 | 乾燥酸化 B/A 前指数 (um/hr) |
| `ox.dry.ae` | 2.00 | 乾燥酸化 B/A 活性化エネルギー (eV) |
| `ox.wet.b0` | 3.86e2 | 湿潤酸化 B 前指数 |
| `ox.wet.be` | 0.78 | 湿潤酸化 B 活性化エネルギー (eV) |
| `ox.wet.a0` | 1.63e8 | 湿潤酸化 B/A 前指数 |
| `ox.wet.ae` | 2.05 | 湿潤酸化 B/A 活性化エネルギー (eV) |

未設定時は元の `arrhenius(7.72e2, 1.23, T)` 等と完全に同じ計算になるため、
このリファクタ単体では出力はビット不変。

### 3.2 圧力・HCl・方位スケーリング

`deal_grove_dry`/`deal_grove_wet`/`deal_grove_step` に
`pressure_atm=1.0, hcl_frac=0.0, orient_factor=1.0` を追加。すべて既定値
でビット不変。

- 圧力: `B *= pressure_atm`, `B/A *= pressure_atm^0.75`
  (Deal & Grove 1965 の高圧酸化拡張)。
- HCl: `hcl_gain = ParamDB.get("ox.hcl.gain", 6.0)` として
  `B *= (1+hcl_gain*hcl_frac)`, `B/A *= (1+hcl_gain*hcl_frac)`。
- 方位: `B/A *= orient_factor` のみ(線形・界面反応律速項のみが方位
  依存、放物線・拡散律速項 B は方位非依存という Deal-Grove の標準的
  取り扱いに従う)。`proc::oxidize()` 側で `orient="<111>"` のとき
  `orient_factor = ParamDB.get("ox.orient.ratio111", 1.68)`、それ以外
  (既定 `"<100>"`)は `1.0`。

### 3.3 Massoud 薄膜増速項

新関数 `deal_grove_step_massoud(x0_um, dt_min, T_c, wet, massoud_c,
massoud_l, pressure_atm, hcl_frac, orient_factor)` を追加。
`massoud_c == 0`(既定)なら `deal_grove_step` に委譲し完全に同じ結果を返す
(ビット不変)。`massoud_c != 0` の場合は成長速度式

```
dx/dt = B/(2x+A) * (1 + C * exp(-x / L))
```

を RK4 で数値積分する(`dt_min` を固定 20 サブステップに分割、時間は B/A
がアレニウス「時/hour」単位で定義されているため hour に変換してから積分)。
`proc::oxidize()` はこの関数を `deal_grove_step` の代わりに呼び、
`ox.massoud.c`(既定 0.0、OFF)・`ox.massoud.l`(既定 0.01 um = 10 nm)を
ParamDB から読む。

**既定 ON/OFF の判断**: 既定 OFF とした。理由は (a) 既存の複数テスト
(`tests/test_oxidize_flow.cpp` の `test_analytic_agreement` など)が
`proc::oxidize()` の既定呼び出し結果と `deal_grove_step()` の厳密一致を
アサートしており、デフォルト物理を変えることは本タスクのスコープ外の
広範な回帰リスクを生む、(b) `docs/IMPLEMENTATION_PLAN_v2.md` C-3 の DoD が
明示的に「厚膜域は現行 DG とビット一致(補正 off)」を要求している。
校正済みの値 (`ox.massoud.c=0.9`, `ox.massoud.l=0.01`) はいつでも
ParamDB 経由(C++/Python/deck の `pdbset`)で有効化できる。

### 3.4 `proc::oxidize()` の新規オプション引数

```cpp
double oxidize(SimState& st, double time_s, double temp_k, bool wet = false,
               std::ostream* log = nullptr, double pressure_atm = 1.0,
               double hcl_frac = 0.0, const std::string& orient = "<100>");
```

既定値はすべて現行(pre-C-3)挙動を再現する。`log`, `pressure_atm`,
`hcl_frac`, `orient` は末尾に追加したため、既存の全呼び出し箇所
(位置引数を末尾まで渡していない)は無変更でコンパイル可能。

Python/pybind バインディングは CLAUDE.md 必須のトリプルに従い追加:
`python/_cprocess.cpp` の `proc_oxidize` に `pressure_atm`/`hcl_frac`/
`orient` の `py::arg`(既定値付き)、`python/cprocess/simulation.py` の
`Simulation.oxidize()` に同名キーワード引数を追加して委譲。

## 4. どう検証するか(結果の目視・抽出手段)

- **エンジニアがこのステップの結果を確認する方法**:
  1. `proc::oxidize()` の戻り値(新オイド厚 cm)と、`log`(`std::ostream*`)
     に出力される `[oxidize] dry|wet ... K ... s: tox X -> Y um` の行
     (P1-6 由来、変更なし)で、厚み推移を追える。
  2. Massoud/圧力/HCl/方位の効果は `deal_grove_step`/`deal_grove_dry`/
     `deal_grove_wet` を直接呼んで(`tests/test_oxidation.cpp` のように)
     `DealGroveParams{A,B,tau}` や成長厚を比較することで定量確認できる。
  3. ParamDB の現在の上書き値一覧は `proc::pdbget_all()`
     (`ParamDB::instance().all()`)で取得可能(deck の `pdbset`/`pdbget`
     コマンドや Python `Simulation` 経由でも到達可能; 既存の P1-10 機構
     をそのまま利用)。
  4. `proc::save`/`export_device` で書き出した VTU をビジュアライザで
     見れば、酸化膜厚の空間分布(1D ブランケットなので厚さは一様)を
     目視確認できる(既存の可視化パス、変更なし)。

## 5. テスト

- `tests/test_oxidation.cpp` に追加:
  - Massoud opt-in で 900℃/10nm 域の成長が pure Deal-Grove の 1.10x 超
    (旧 `test_sprocess_parity.cpp` の C-3 チェックを移設)。
  - 既定(ParamDB 未設定)での `proc::oxidize()` と `deal_grove_step()`
    のビット一致(リファクタが無害であることの回帰アンカー)。
  - 圧力スケーリング: 5atm > 1atm の成長厚。
  - HCl スケーリング: 3% HCl > HCl なしの成長厚。
  - 方位スケーリング: `<111>` の B/A が `<100>` の 1.68x(±5%)、および
    実際の成長厚が `<111>` > `<100>`。
- `tests/test_benchmarks.cpp` の Tier A 酸化ベンチマーク 5 件
  (dry 1000℃/120min, dry 1100℃/30・60min, wet 1000℃/30・60min)は
  無変更で PASS を維持(Massoud 既定 OFF のため、これらは元の
  `deal_grove_step` 経路と完全に同じ計算)。
- `python/test_comprehensive.py` に `oxidize(pressure_atm=..., hcl_frac=...,
  orient=...)` のスモークテストを追加。

## 6. やらないこと

- Massoud 既定値を ON に変更すること(スコープ外、上記 3.3 参照)。
- `oxidize_2d`(LOCOS)への圧力/HCl/方位/Massoud の配線(ブランケット
  `oxidize()` のみが対象。2D 版は別タスクとする)。
- HCl による塩素の Si 中への取り込みや金属イオンゲッタリングなど、
  成長速度以外の HCl 酸化の副次効果のモデル化。
- 圧力容器/LPCVD 装置固有の伝熱・ガス流モデル(本タスクは Deal-Grove
  係数への一様スケーリングのみ)。
