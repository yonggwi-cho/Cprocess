# P2-8: ドーパント追加 (In, C, F, Ge)

## 目的

現行の 4 種 (B/P/As/Sb) に In・C・F・Ge を追加する。In は低速・
低固溶度のアクセプタ (チャネルプロファイル制御用)、C は C-I クラスタで
TED を抑制する中性種 (SiGe:C 技術の中核)、F は中性の高速拡散種、
Ge は不動のマーカー/歪み種 (SiGe の入口)。**中性種の導入**が本タスクの
構造的変更点: 電荷中性計算とフィールド増速から除外する仕組みを入れる。

## 現状コード

- `src/materials.cpp` — `kDopants` テーブル (name, symbol, DopType,
  z, m, d0/e0, dm/em, dmm/emm, dp/ep, ss_pre/ss_e, fi, range table)
- `include/cprocess/materials.hpp` — `enum class DopType { donor, acceptor }`
- `src/diffusion.cpp` — DopType を参照する箇所:
  (1) `run()`/`run_ted()` の nni ループ
  `nnet += (type == donor) ? c : -c`、(2) field_enh 分岐
  `(donor && ntype) || (acceptor && !ntype)`
- `src/mc_implant.cpp` — ドーパントは z/m のみ使用 (追加変更不要のはず —
  確認して仕様に記録)
- P2-1 (フル点欠陥) がマージ済みなら C-I シンクは P2-1 の反応項に、
  未マージなら現行 `run_ted` に入れる — **本仕様は現行 run_ted 前提**で
  書き、P2-1 側の移植規則を「やらないこと」の注記に残す

## 実装手順

1. `DopType` に `neutral` を追加 (`materials.hpp`)。
   `diffusion.cpp` の 2 箇所を修正:
   - nni ループ: `if (type == neutral) continue;`
   - field_enh: neutral は増速なし (既存分岐に該当しないので
     変更不要だが、明示の `if (dp.type == DopType::neutral)` スキップを
     入れて意図を明確化)
2. `kDopants` に 4 行追加 (教科書レベル・オーダー重視の既定値):
   ```
   In: acceptor, Z=49, m=114.82
       d0=0.785, e0=3.63 (Fair)、dp=0.415, ep=3.63
       ss: 5e17 @1000C 相当 → ss_pre=6.9e20, ss_e=0.78
       fi=0.2 (空孔優勢)
       range: Sb に近い {10:9nm,4nm}..{200:100nm,34nm} (Sb 表を 1.05 倍)
   C:  neutral, Z=6, m=12.011
       d0=0.95, e0=3.04 (置換 C)
       ss_pre=4e24, ss_e=1.0 (低固溶度 ~9e19@1000C)
       fi=1.0 (I 媒介)
       range: B 表の深さ 0.8 倍
   F:  neutral, Z=9, m=18.998
       d0=1.0e-2, e0=2.2 (高速; 実際は欠陥依存だが定数近似)
       ss_pre=1e23, ss_e=0.8, fi=0.5
       range: B と P の中間 (B 表 0.7 倍)
   Ge: neutral, Z=32, m=72.63
       d0=2.5e3? → 不動扱い: d0=0.0 (D=0; 拡散ソルバーは D=0 種を
       スキップ — dcell=0 のセルは面フラックスが立たないので自然に不動。
       全セル D=0 の種のソルブが退化しないことを確認し、退化するなら
       種スキップ条件 `if (max dcell == 0) continue;` を run() に追加)
       ss_pre=0 (クランプなし), fi=0
       range: As 表と同一
   ```
   数値は仕様値として固定 (較正は ParamDB (P1-10) で上書き可能)
3. C-I シンク (TED 抑制) — `run_ted()` に追加:
   - C 場が存在する場合、ψ (I 過剰) の更新後に
     `dψ = −k_CI · C_C · ψ · dt` (陽的、サブサイクル不要な安定条件
     `k_CI·C_C·dt < 0.5` を dt クランプで保証) を適用し、
     同量を不動場 `"C_cl"` に積算 (C 自体は減らさない — C-I ペアの
     C は再放出される近似; 仕様として明記)
   - `k_CI` は ParamDB キー `"ted.k_ci"` 既定 2e-21 cm³/s
     (C 1e19 cm⁻³ で時定数 ~50 s になるオーダー)
4. MC 注入: `find_dopant` が返る種はそのまま动く (z/m 使用のみ) —
   4 種で apply_mc_implant がエラーなく走ることをテストで確認
5. Python/pybind: 変更不要 (種名文字列で流れる)。テストのみ追加

## テスト仕様

`tests/test_new_dopants.cpp` 新設 (CMakeLists foreach 登録):

1. **基本動作**: 各新種について implant_gauss (1e14, rp=50nm) →
   diffuse (10min/1000C) がエラーなく完走、dose 変化 < 0.5%、
   field finite
2. **MC 動作**: In と Ge の MC 注入 (30keV, 2e4 ions) が deposited > 0
3. **中性種の電気的中立**: P 1e18 一様 + Ge 1e20 背景ありなしで
   P の拡散結果が一致 (相対差 < 1e-9) — nni に Ge が入っていない証明
4. **Ge 不動**: Ge プロファイルが diffuse 前後で不変 (< 1e-12)
5. **C の TED 抑制**: B 1e14 (damage=true) + C 1e19 背景の run_ted と
   C なしの run_ted を比較 → C ありの B 拡がり増分が 20% 以上小さい。
   `"C_cl"` 場が生成され積算量 > 0
6. **In の低速性**: 同条件で In の拡がり < B の拡がり (1000C 30min)
7. 既存テスト全 PASS (B/P/As/Sb の挙動不変)

Python (`python/test_comprehensive.py`): 4 種の implant+diffuse スモーク
(1 テスト関数にまとめてよい) + C 抑制の定性チェック。

## 完了条件 (DoD)

- [ ] DopType::neutral 追加、nni/field_enh の除外 2 箇所
- [ ] kDopants 4 行 (In/C/F/Ge) + range 表
- [ ] C-I シンク (run_ted) + "C_cl" 積算 + ParamDB キー
- [ ] 上記 7 C++ テスト + Python テスト PASS
- [ ] コミットメッセージに `P2-8` を含める

## やらないこと

- P2-1 (フル I/V) への C-I 反応の移植 — P2-1 実装時に
  「C-I シンクを I 方程式の反応項へ移す」1 行を P2-1 側 DoD に追加する
- Ge 濃度→格子歪み (P3 SiGe タスク)
- F の欠陥依存拡散 (vacancy trap モデル)
- In-C 共注入等の相互作用
