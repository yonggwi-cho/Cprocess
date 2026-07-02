# S-5: 適応時間刻み (step-doubling 局所誤差制御)

## 目的

拡散/TED の時間刻みは現在 `time/50` (または opts.dt) の等間隔で、TED 初期の
急峻な過渡には過大、後半には過小になる。step-doubling (1 回の dt ステップと
2 回の dt/2 ステップの比較) による局所誤差推定で dt を自動制御する。

**コスト注意**: 採択 1 ステップあたり線形ソルブが約 3 倍になる (dt 1 回 +
dt/2 2 回)。過渡で大きく刻み幅を伸ばせる問題 (TED、RTA スパイク) でのみ
得をする設計であり、デフォルトは **無効** (完全後方互換)。

内部ソルバー変更のため Python 表面は不要 (DiffuseOpts の新フィールドは
デフォルト値で既存呼出しに影響しない)。

## 現状コード

- `src/diffusion.cpp` — `DiffusionSolver::run` / `run_ted`:
  `dt0 = (o.dt > 0) ? min(o.dt, o.time) : o.time / 50.0;`
  `nsteps = ceil(time/dt0)` の固定ループ。各ステップは
  `cold[s] = *fields[s].conc` を保存 → Picard 反復 → 負値クランプ、の構造。
  run_ted はさらにステップ冒頭で psi を 1 回線形陰解 (`psi_old = psi` 保存)
- `include/cprocess/diffusion.hpp` — `DiffuseOpts { double dt = 0; ... }`

## 実装手順

1. **DiffuseOpts 拡張** (diffusion.hpp):
   ```cpp
   bool adaptive_dt = false;   // step-doubling dt control (only when dt == 0)
   double dt_tol = 0.05;       // relative local-error tolerance
   // Optional instrumentation: if non-null, accepted dt values are appended.
   mutable std::vector<double>* step_log = nullptr;
   ```
   **有効条件を明確に**: 適応刻みは `o.dt == 0 && o.adaptive_dt == true` の
   ときだけ。`o.dt > 0` または `adaptive_dt == false` (デフォルト) なら
   現行の固定刻み経路を**一切変更せず**通す (既存結果 bit 一致)。
2. **1 ステップの試行を関数化**: run / run_ted のステップ本体 (cold 保存の
   直後〜負値クランプまで。run_ted では psi 陰解 + S 計算 + dopant Picard の
   全体) を private ヘルパーに切り出す:
   ```cpp
   // Advances all fields (and psi for the TED variant) by one backward-Euler
   // step of length dt, starting from the state currently held in the field
   // vectors. Overwrites the fields with the end-of-step state.
   void step_once(...);      // run 用
   void step_once_ted(...);  // run_ted 用 (psi を含む)
   ```
   引数はステップ本体が現在参照しているローカル (cold, dcell, nni, rhs, x,
   grad, bcface, dt, opts, ni, cfloor, TED では psi/psi_bcface/dI/cstar/
   k_rec) を素直に渡す。**中身の演算列は移動のみで変更しない** (固定刻み
   経路の bit 一致を守るため)。
3. **適応ループ** (run / run_ted 双方の `for (step ...)` を置換する分岐):
   ```
   dt = time / 50            // 初期刻み
   dt_min = time / 10000
   t = 0
   while t < time - 1e-12 * time:
     dt = min(dt, time - t)
     saved = 全フィールドのコピー (TED では psi も)
     // 試行 1: dt を 1 回
     step_once(dt)           → c_dt に退避、フィールドを saved に戻す
     // 試行 2: dt/2 を 2 回
     step_once(dt/2); step_once(dt/2)   → フィールドが c_half を保持
     err = max over 全種・全セル |c_dt - c_half| / (|c_half| + cfloor)
           (cfloor = 1e-3 * max(全種の max c_half, 1.0);
            TED では psi も同じ式で err に含める)
     if err > o.dt_tol:
       // 棄却: フィールドを saved に戻し dt /= 2
       if dt/2 < dt_min: throw std::runtime_error(
           "diffusion: adaptive dt underflow (dt < time/10000)")
       dt /= 2; continue
     // 採択: c_half (2 半ステップの解、より高精度) を採用
     t += dt
     if o.step_log: o.step_log->push_back(dt)
     if err < o.dt_tol / 4: dt = min(dt * 1.5, time / 10)
     dt = min(dt, time - t)
   ```
   - 採択解は **c_half 側** (フィールドはそのまま)。c_dt は誤差推定にのみ使う
   - 棄却時のフィールド復元は saved からの単純代入
   - ログ (`verbosity>=1`) は「採択ステップごと」に現行フォーマットで
     `t` と dt を出す (10 分の 1 間引きは `nsteps` が不定なので廃止し、
     採択ステップ数で間引くか毎回出力でよい — 適応経路のみの話)
4. **質量・peak の最終ログ**は現行コードを共用 (ループ後は同一)。
5. **P1-5 (RTA ランプ) への前方互換ノート** (仕様書コメントとして): ランプ
   温度テーブル導入時は「コントローラが提案した dt を次のブレークポイント
   までにクランプ」する 1 行を採択前に挟むだけでよい。本タスクでは実装しない。

## テスト仕様

新規 `tests/test_adaptive_dt.cpp`、CMakeLists の foreach に `adaptive_dt` を
追加。セットアップは tests/test_ted.cpp の TED ケース (Gaussian B +
"+1" psi シード、`make_box_mesh`) を踏襲。

1. **過渡での刻み挙動**: TED ケース (急峻な初期過渡) を
   `adaptive_dt=true, dt=0, dt_tol=0.05, step_log=&log` で `run_ted`。
   assert: `log.front() < log.back()` (初期の採択 dt < 最終の採択 dt)、
   かつ `log.size() >= 3`。
2. **精度と効率**: 同ケースを固定刻み `dt = time/500` で解いた参照と比較。
   assert: (a) 各セルの最終 B 濃度の相対差
   `max |c_ad - c_ref| / (peak_ref)` **< 2%** (2e-2)、
   (b) 適応の採択ステップ数 `log.size()` **< 250** (参照の 500 の半分未満)。
3. **後方互換 bit 一致**: 同ケースを `adaptive_dt=false` (デフォルト) で
   実行した結果が、本タスク適用前のコードと同一経路 (コードレビューで確認)
   であることに加え、`opts.dt = time/50` 明示指定の実行と全セル `==`
   (bit 一致) すること。step_once への切出しが演算列を変えていないことの
   検証になる。
4. **dt アンダーフロー**: `dt_tol = 1e-12` (到達不能) で
   `std::runtime_error` が投げられること。
5. **run 側**: 通常の B アニール (TED なし) でもテスト 2 と同基準
   (相対差 < 2%、ステップ数 < 250) で `run` の適応経路を確認。

## 完了条件 (DoD)

- [ ] `DiffuseOpts::adaptive_dt / dt_tol / step_log` 追加 (デフォルト無効)
- [ ] run / run_ted 双方に step-doubling 制御が入り、固定刻み経路は bit 不変
- [ ] `tests/test_adaptive_dt.cpp` の 5 項目 PASS、CMakeLists 登録
- [ ] 既存全テスト PASS (デフォルト無効なので変化しないはず)
- [ ] コミットメッセージに `S-5` を含める

## やらないこと

- 埋込み Runge-Kutta ペア・二次外挿など step-doubling 以外の誤差推定
- 適応刻みのデフォルト有効化 (較正が済むまで opt-in)
- RTA 温度ランプ (P1-5) 本体 — 前方互換ノートのみ
- Picard 反復数に基づく刻み制御 (誤差ベースのみ)
- Python API への露出 (proc:: / Simulation への opts 追加は P1-5 側で一括)
