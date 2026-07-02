# P2-2: ドーパントクラスタ (BIC / As-V) — 活性化の時間依存

## 目的

現行の活性化は固溶度クランプ (P1-3) のみで瞬時・可逆であり、逆アニール
(低温で活性率が一旦**下がる**) や高濃度 As の不活性化を再現できない。
本タスクで不動・電気的不活性の**クラスタ場**を導入する:

- B: BIC (B_m I_n クラスタ; 有効組成 m=3, n=1 の単一場で近似) — 過剰 I が
  あるほど形成が速い → TED 中の B 不活性化 (逆アニール) を再現
- As: As₄V — 高濃度 As の不活性化を再現

**P2-1 (フル点欠陥) と P1-3 (活性化)・P1-10 (ParamDB) の完了が前提。**

## 現状コード

- `src/diffusion.cpp` — `run_ted` (P2-1 完了後: C_I/C_V/C_311 分割解法、
  反応サブサイクルの骨組みがある)。電荷中性ループ (`nni` 計算) は
  総濃度をそのまま使っている。P1-3 完了後は活性濃度
  `C_act = min(C, C_ss(T))` を電荷中性に使う機構がある。
- `include/cprocess/materials.hpp` — `Dopant` (symbol "B"/"As"、`fi`)、
  `solid_solubility`。
- `st.fields` — `std::map<std::string, std::vector<double>>`。新フィールド名は
  proc 層が自由に追加できる (save() は全フィールドを VTU に出す)。
- P1-5 (RTA ランプ) 完了済みなら 2 温度ステップはランプで書ける。未完なら
  `diffuse_ted` を 2 回呼ぶ。

## 実装手順

1. **クラスタ場**: `st.fields["B_cl"]` / `st.fields["As_cl"]` (cm⁻³、
   クラスタに取り込まれた**ドーパント原子**の濃度)。不動 (拡散させない)、
   電気的不活性。`proc::diffuse_ted` がドーパント場を集める際、シンボル末尾
   `_cl` の場は SpeciesField にしない (拡散対象から除外)。
2. **簡略化した交換動力学** (反応サブサイクル内、C_I/C_V/C_311 と同じ
   前進オイラー・サブステップで一緒に進める):
   - **B (BIC)**:
     ```
     dC_Bcl/dt = k_f_B · C_B_act² / C_ref · (C_I / C_I*) − k_r_B · C_Bcl
     dC_B/dt   = −dC_Bcl/dt        (C_B は可動 B の総濃度)
     ```
     forward は p=2 (実効 B₃I の核形成が B 対から始まる近似)、I 過剰係数
     q=1。`C_ref = 1e20 cm⁻³` (次元合わせの規格化濃度、固定)。
     `k_f_B = params.get("cl.b.kf", 1e-3)` [1/s]、
     `k_r_B = params.get("cl.b.nu0", 1e13) · exp(−params.get("cl.b.eb", 3.6)/kT)`
     [1/s] (Arrhenius 解離; E_b≈3.6 eV → 700 °C で ~1e-6/s ≒ 不解離、
     900 °C で ~1e-2/s ≒ 数分で解離)。
   - **As (As₄V)**: 真の C_As⁴ 律速は極端に硬い (stiff) ため**実効 2 次**で
     近似する (コメントに「As₄V の 4 次核形成を 2 次で実効近似。真の 4 次は
     サブサイクル刻みを桁で縮めるため採らない」と明記):
     ```
     dC_Ascl/dt = k_f_As · C_As_act² / C_ref · (C_V / C_V*) − k_r_As · C_Ascl
     ```
     `cl.as.kf` 既定 5e-4 [1/s]、`cl.as.nu0` 1e13、`cl.as.eb` 3.2 eV。
   - 高濃度閾値: forward 項は `C_act > 0.1 * solid_solubility(d, T)` の
     セルのみで評価 (低濃度での偽クラスタリング防止; テスト 3 の根拠)。
   - サブサイクル刻み `dt_react` の rate_max に `k_f·C_act/C_ref·(C_I/C_I*)`
     と `k_r` を加える (P2-1 の式に追加)。
   - クラスタ形成/解離は点欠陥も交換する (B₃I は I を 1/3 個/B 原子):
     `dC_I/dt += −(1/3)·dC_Bcl/dt`、As₄V は `dC_V/dt += −(1/4)·dC_Ascl/dt`。
3. **活性濃度**: `C_act = max(C_total − C_cl, 0)` をまず計算し、その上で
   P1-3 の固溶度クランプを適用 (`C_act = min(C_act, C_ss)`)。この C_act を
   (a) 電荷中性 (nni ループ)、(b) 上記 forward 項、(c) save() の `*_active`
   出力、に使う。**拡散するのは可動場 C_B (= total − cluster) のまま**;
   実装では `st.fields["B"]` を可動+活性の総和とし、クラスタ分は "B_cl" に
   移し替える (B + B_cl = 全 B)。
4. **結線**: P2-1 の反応サブサイクル関数にクラスタ項を追加するだけで、
   新しい proc:: 関数は不要 (`diffuse_ted` が自動で扱う)。B または As の
   場が存在するとき対応する `*_cl` 場を lazily 作成する。
   pybind/Simulation の追加も不要だが、`sim.field("B_cl")` の取得を
   Python テストで確認する。
5. デッキ/ログ: `[ted]` サマリ行に `B_cl_frac=` (クラスタ率) を追加。

## テスト仕様

`tests/test_dopant_clusters.cpp` (新規、CMakeLists foreach に
`dopant_clusters` 追加)。メッシュは test_ted.cpp の setup を流用。

1. **逆アニールの兆候**: B を高濃度 (dose 1e15, Rp 0.05 µm → ピーク
   ~2e20 cm⁻³) で damage=true 注入。
   - Step 1: `diffuse_ted(10 s, 700 °C)` → 活性率
     `sum(C_act·V)/sum((C_B+C_Bcl)·V)` が **< 0.9** (クラスタ形成で低下)
   - Step 2: 続けて `diffuse_ted(10 min, 900 °C)` (P1-5 のランプがあれば
     ランプ 1 回でもよい) → 活性率が **> 0.95** に回復 (解離)
2. **B 総量保存**: テスト 1 の全工程で `sum((C_B + C_Bcl)·V)` の変化が
   相対差 < 0.1% (交換は厳密に対で行われる)。
3. **低濃度で偽クラスタなし**: dose 1e12 (ピーク ~2e17) の B を damage=true
   注入 → `diffuse_ted(60 s, 800 °C)` → クラスタ率 < 1%。
4. **As₄V**: As を dose 1e16 (ピーク > 1e21) で注入、`diffuse_ted(60 s,
   900 °C)` → As クラスタ率 > 5%、As 総量 (As + As_cl) 保存 < 0.1%。
5. **既存テスト**: test_ted / test_point_defects (P2-1) が PASS すること
   (基準の増速テストは dose 1e14 ≒ ピーク 2e19 でクラスタ閾値近傍のため、
   もし 1.3x を割る場合は `cl.b.kf` 既定値を下げて較正してよい — その場合
   テスト 1 の基準も同時に満たすこと)。

Python (`python/test_comprehensive.py`): 高濃度 B damage 注入 →
`diffuse(ted=True)` → `sim.field("B_cl").sum() > 0`、
`(sim.field("B") + sim.field("B_cl"))` の dose が注入 dose と < 1% で一致。

## 完了条件 (DoD)

- [ ] "B_cl"/"As_cl" 場と交換動力学 (上記の具体式・ParamDB キー) 実装
- [ ] 電荷中性・save の active 出力が `total − cluster` ベース
- [ ] C++ テスト 5 本 PASS、CMakeLists 登録
- [ ] Python テスト追加、`python3 python/test_comprehensive.py` PASS
- [ ] 総量保存 < 0.1% (テスト 2/4)
- [ ] コミットメッセージに `P2-2` を含める

## やらないこと

- クラスタサイズ分布 (B₂I, B₃I₂ … の個別追跡) — 単一有効場のみ
- 真の 4 次 As 動力学 (実効 2 次で近似、コメント必須)
- P/Sb のクラスタ (必要になったら同じ機構にテーブル追加)
- アモルファス化領域の特別扱い / SPER (P3)
- クラスタの陰解法・Newton (S-3 後に再検討)
