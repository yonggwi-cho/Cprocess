# P1-2: MC 損傷 → TED 接続 (Kinchin-Pease 損傷場で格子間原子をシード)

## 目的

TED の格子間原子シードは現在 "+1" モデル (`seed_interstitials`,
src/process.cpp): 注入で増えたドーパント濃度分布をそのまま I 場へコピーする。
一方 MC 注入 (src/mc_implant.cpp) は Kinchin-Pease 変位密度 [cm^-3] を
物理的に計算しているのに、`proc::implant_mc` は damage 出力を受け取らず捨てて
いる (`apply_mc_implant` の `damage_conc` 引数に nullptr / デフォルト)。
MC 注入では MC 自身の損傷場を Frenkel 対生存率でスケールして I 場をシード
するよう接続する。損傷はドーパント静止点より浅い側にピークを持つため、
"+1" より物理的に正しい深さ分布になる。解析 Gaussian 注入は損傷情報を
持たないので "+1" のまま。

## 現状コード

- `src/mc_implant.cpp` — `apply_mc_implant(mesh, silicon_mask, p, conc,
  damage_conc)`:
  - `damage_conc` 非 null **かつ `p.channeling == true`** のときだけ、末尾で
    `(*damage_conc)[i] = dmg_store.counts[i] * weight / mesh.cell_vol[i]`
    (変位原子密度 [cm^-3]) を書き込む。channeling=false では DamageStore 自体
    が作られず damage は得られない
  - Kinchin-Pease: 反跳エネルギー t_recoil >= kEdSi (15 eV) で
    ndis = max(1, t_recoil/(2*kEdSi)) を蓄積
  - アモルファス化閾値 `kNamorph = 6.25e21` [cm^-3] が anonymous namespace の
    ローカル定数として定義され、`f_amor = dens / kNamorph` に使われている
- `src/process.cpp` — `proc::implant_mc`:
  - 非スタック経路: `apply_mc_implant(st.mesh, silicon_mask(st), p, f)`
    (damage_conc 省略 = nullptr)
  - スタック経路 (st.has_stack): `apply_mc_implant(st.stack, stack_si, p,
    stack_conc, nullptr)` → `transfer_field_nearest(st.stack, stack_conc,
    st.mesh)` で作業メッシュへ転写
  - どちらも `seed_damage` なら `seed_interstitials(st, before, f)` ("+1")
- `src/process.cpp` — `seed_interstitials(st, before, after)`:
  I[i] += max(0, after[i]-before[i])
- `include/cprocess/field_transfer.hpp` — `transfer_field_nearest(src_mesh,
  src_field, dst_mesh)` (スタック→作業メッシュの転写に既用)
- `src/diffusion.cpp` — `run_ted` は st.fields["I"] (= psi) を消費。変更不要
- pybind `proc_implant_mc` の既定は `channeling=false`、
  `Simulation.implant` の既定も `channeling=False`

## 実装手順

1. `include/cprocess/materials.hpp` に定数を 2 つ追加 (namespace cp 直下):
   ```cpp
   // Fraction of Kinchin-Pease Frenkel pairs surviving in-cascade
   // recombination; the survivors seed the excess-interstitial field for TED
   // (cf. "+1" model: net excess ~ dose, i.e. ~1% of total displacements).
   constexpr double kFrenkelSurvival = 0.01;

   // Displacement density at which crystalline Si is fully amorphized
   // [cm^-3] (~12.5% of the atomic density). Also caps the seeded
   // interstitial excess: cells at the cap are amorphous and their TED
   // physics differs, but capping keeps the free supersaturation bounded.
   constexpr double kAmorphizationDensity = 6.25e21;
   ```
2. `src/mc_implant.cpp` — ローカル定数 `kNamorph` を削除し、参照 2 箇所
   (`f_amor` の割り算、コメント) を `kAmorphizationDensity` に置き換える
   (`materials.hpp` は mc_implant.hpp 経由で既に include されている)
3. `src/process.cpp` — anonymous namespace に新ヘルパーを追加:
   ```cpp
   // Seed the interstitial excess from the MC Kinchin-Pease damage field:
   //   I += kFrenkelSurvival * damage, capped at kAmorphizationDensity.
   void seed_interstitials_from_damage(SimState& st,
                                       const std::vector<double>& damage) {
     auto& I = st.fields["I"];
     I.resize(st.mesh.cells.size(), 0.0);
     for (std::size_t i = 0; i < damage.size() && i < I.size(); ++i) {
       I[i] += kFrenkelSurvival * damage[i];
       I[i] = std::min(I[i], kAmorphizationDensity);
     }
   }
   ```
4. `proc::implant_mc` の配管を変更する。方針: `seed_damage && p.channeling`
   のときは損傷ベースでシードし、`seed_damage && !p.channeling` のときは
   従来の "+1" にフォールバック (channeling off では MC 損傷が存在しないため)。
   - **非スタック経路**:
     ```cpp
     std::vector<double> dmg;
     const bool use_damage = seed_damage && p.channeling;
     const McImplantStats s = apply_mc_implant(
         st.mesh, silicon_mask(st), p, f, use_damage ? &dmg : nullptr);
     if (use_damage)       seed_interstitials_from_damage(st, dmg);
     else if (seed_damage) seed_interstitials(st, before, f);
     ```
     (`before` の取得は従来どおり seed_damage 時のみ)
   - **スタック経路**: 同様に `stack_dmg` を確保して
     `apply_mc_implant(st.stack, stack_si, p, stack_conc,
     use_damage ? &stack_dmg : nullptr)` とし、use_damage なら
     `const std::vector<double> dmg_transferred =
     transfer_field_nearest(st.stack, stack_dmg, st.mesh);` →
     `seed_interstitials_from_damage(st, dmg_transferred)`。
     それ以外は従来の "+1"
   - ログ: use_damage のとき
     `[implant] damage seed: peak I=... cm^-3 (KP damage x 0.01)`、
     フォールバック時 `[implant] damage: '+1' model (channeling off)` を 1 行
5. `implant_gauss` は変更しない ("+1" のまま)。`seed_interstitials` は残す
6. ヘッダ (`include/cprocess/process.hpp`) の `implant_mc` コメントを更新:
   「seed_damage=true かつ channeling=true なら Kinchin-Pease 損傷 ×
   kFrenkelSurvival (上限 kAmorphizationDensity) を I 場へシード。
   channeling=false は "+1" フォールバック」
7. proc:: 関数のシグネチャは不変なので pybind の変更は不要だが、
   `python/_cprocess.cpp` の `proc_implant_mc` docstring と
   `python/cprocess/simulation.py` の `Simulation.implant` docstring に
   「mc=True かつ channeling=True かつ damage=True で MC 損傷ベースの
   I シード (それ以外の damage=True は "+1")」を追記する。
   新しい proc:: 関数は追加しないため新規バインディングは不要 (CLAUDE.md
   の規約はシグネチャ不変のため docstring 更新 + Python テストで満たす)

## テスト仕様

`tests/test_ted.cpp` に追加 (既存 3 テストの後。質量加重の平均深さヘルパー
`mean_depth(st, c) = Σ c_i V_i (ztop - z_i) / Σ c_i V_i` を追加する):

1. **損傷ピークは Rp より浅い**: `proc::mesh_box(st, 0, 0.3e-4, 0, 0.3e-4,
   0, 0.5e-4, 4, 4, 50)`、B 50 keV, dose=1e14, ions=200000, threads=1,
   channeling=true, seed_damage=true で `proc::implant_mc`。
   `CHECK(st.fields.count("I") == 1)` かつ
   **`mean_depth(st, I) < mean_depth(st, B)`** を CHECK
   (核阻止能は静止点より手前で最大 → 損傷は浅い)
2. **I 積分はドーズに線形**: 同条件で dose=1e14 と dose=2e14 (同 seed,
   ions=200000, threads=1) をそれぞれ新しい SimState で実行し、
   I の総量 (Σ I_i V_i) の比が **[1.8, 2.2]** (2 倍 ±10%) に入ることを CHECK
   (channeling の損傷フィードバックで完全線形にはならない)
3. **MC シードでも TED が効く**: 既存テスト 1 と同じ比較を MC 版で行う。
   equilibrium: B 20 keV MC 注入 (dose=1e14, ions=100000, threads=1,
   channeling=true, seed_damage=false) → `proc::diffuse`
   (temp=1173.15, time=60)。TED: 同一注入を seed_damage=true で行い
   `proc::diffuse_ted`。**spread_ted > 1.3 * spread_eq** を CHECK
4. **シード上限**: 極端条件 dose=1e16, B 50 keV, ions=100000, threads=1,
   channeling=true, seed_damage=true で
   **`max(I) <= kAmorphizationDensity * (1 + 1e-12)`** を CHECK
5. **既存挙動の非退行**: 既存の test_ted 3 テスト (Gaussian "+1" ベース) が
   無変更で PASS すること。また channeling=false + seed_damage=true の MC
   注入で I 場が生成される ("+1" フォールバック、Ipeak > 0) ことを CHECK

`python/test_comprehensive.py` に追加:

- `r = sim.implant("B", dose=1e14, energy=50, mc=True, ions=50000,
  channeling=True, damage=True, threads=1, seed=1)` 後、
  `"I" in sim.field_names()`、`sim.field("I").max() > 0`、
  `sim.field("I").max() <= 6.25e21` を assert
- I の平均深さ < B の平均深さ (セル体積 × 濃度で加重、
  `sim.cell_centroids[:,2]` 使用) を assert
- その後 `sim.diffuse(time=1, temp=900, ted=True)` が完走することを確認

実行: `cmake --build build -j$(nproc) && ctest --test-dir build` 全 PASS、
`build/_cprocess*.so` を `python/cprocess/` へコピーして
`python3 python/test_comprehensive.py`。

## 完了条件 (DoD)

- [ ] `kFrenkelSurvival` / `kAmorphizationDensity` が materials.hpp の
      doc コメント付き定数になり、mc_implant.cpp の `kNamorph` が
      `kAmorphizationDensity` に一本化されている
- [ ] `proc::implant_mc` が非スタック/スタック**両経路**で damage 配列を
      受け取り、スタック経路は `transfer_field_nearest` で転写してからシード
- [ ] seed_damage && channeling → 損傷シード (×0.01, cap 付き)、
      seed_damage && !channeling → "+1" フォールバック (ログで区別可能)
- [ ] Gaussian 注入 ("+1") と既存 test_ted テストは無変更で PASS
- [ ] test_ted.cpp に上記 5 ケース追加、Python テスト追加、全テスト PASS
- [ ] コミットメッセージに `P1-2` を含める

## やらないこと

- I と V の 2 種輸送、{311} クラスタの明示的モデル — P2-1。
  run_ted の物理 (kSmax 緩衝、表面シンク) は一切変更しない
- アモルファス層の別扱い (SPER、アモルファス中の拡散) — cap のみで近似
- kFrenkelSurvival の温度/ドーズ依存化、較正 — 定数 0.01 固定
  (上書き機構は P1-10)
- Gaussian 注入への損傷モデル追加 (損傷情報が存在しない)
- channeling=false の MC に損傷蓄積を追加すること (決定性を壊すため現状維持)
- 新しい proc:: 関数・デッキコマンドの追加 (既存シグネチャのまま)
