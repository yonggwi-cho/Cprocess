"""
ボロン注入 + 熱拡散シミュレーション
=====================================
Si(100) 基板に B を 50 keV / 1e13 cm^-2 で注入し、
1000°C 30分のアニールを行う。
"""

import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'python'))

import numpy as np
import cprocess as cp

# --------------------------------------------------------------------------
# 1. メッシュ作成  (0.5 x 0.5 x 1.2 µm, z 方向が深さ)
# --------------------------------------------------------------------------
um = cp.um   # 1e-4 cm

print("=== メッシュ作成 ===")
st = cp.SimState()
log = cp.run_deck(f"""
mesh box xmax=0.5um ymax=0.5um zmax=1.2um nx=8 ny=8 nz=80

# 背景 P ドープ (1e15 cm^-3)
init species=P conc=1e15
""", st)
print(log)

mesh = st.mesh
nc   = mesh.n_cells
z_cm = mesh.cell_cent[:, 2]        # セル重心 z 座標 [cm]
vol  = mesh.cell_vol                # セル体積 [cm^3]
z_top = z_cm.max()                  # 表面 z 座標

print(f"セル数: {nc},  ノード数: {mesh.n_nodes}")
print(f"メッシュサイズ: {z_top*1e4:.1f} µm (深さ)")

# --------------------------------------------------------------------------
# 2. Monte Carlo イオン注入 (B 50 keV, 1e13 cm^-2)
# --------------------------------------------------------------------------
print("\n=== MC イオン注入 (B, 50 keV, 1e13 cm^-2) ===")
log = cp.run_deck("""
implant species=B energy=50keV dose=1e13 method=mc ions=100000 seed=7
""", st)
print(log)

B_after_implant = st.get_field('B').copy()

# Rp / dRp の確認
depth = z_top - z_cm          # 深さ [cm]
total_atoms = np.sum(B_after_implant * vol)
xy_area = (0.5 * um)**2       # 注入面積 [cm^2]
dose_num = total_atoms / xy_area
Rp_num  = np.sum(B_after_implant * vol * depth) / total_atoms
dRp_num = np.sqrt(np.sum(B_after_implant * vol * (depth - Rp_num)**2) / total_atoms)
print(f"数値 Rp  = {Rp_num*1e7:.1f} nm")
print(f"数値 dRp = {dRp_num*1e7:.1f} nm")
print(f"注入ドーズ: {dose_num:.3e} cm^-2 (理論値 1e13)")

# --------------------------------------------------------------------------
# 3. 熱拡散  (1000°C, 30 分)
# --------------------------------------------------------------------------
print("\n=== 熱拡散 (1000°C, 30 min) ===")
log = cp.run_deck("""
diffuse time=30min temp=1000C dt=30s
""", st)
print(log)

B_after_diffuse = st.get_field('B').copy()

# --------------------------------------------------------------------------
# 4. 結果の集約 (1D 深さプロファイル)
# --------------------------------------------------------------------------
print("\n=== 深さプロファイル (深さ方向平均) ===")

# z の一意な値（各スライス）を取得してプロファイルを作る
z_vals = np.unique(np.round(z_cm, 12))[::-1]   # 深さ降順 (表面→深部)
depths_nm, prof_impl, prof_diff = [], [], []

for zv in z_vals:
    idx = np.where(np.abs(z_cm - zv) < 1e-8)[0]
    if len(idx) == 0:
        continue
    w = vol[idx] / vol[idx].sum()
    depths_nm.append((z_top - zv) * 1e7)         # nm
    prof_impl.append(np.sum(B_after_implant[idx] * w))
    prof_diff.append(np.sum(B_after_diffuse[idx] * w))

depths_nm = np.array(depths_nm)
prof_impl  = np.array(prof_impl)
prof_diff  = np.array(prof_diff)

# テキスト出力
print(f"\n{'深さ[nm]':>10}  {'注入直後[cm^-3]':>18}  {'拡散後[cm^-3]':>18}")
print("-" * 52)
step = max(1, len(depths_nm) // 20)
for i in range(0, len(depths_nm), step):
    print(f"{depths_nm[i]:10.1f}  {prof_impl[i]:18.3e}  {prof_diff[i]:18.3e}")

# ピーク / 接合深さ
P_bg = 1e15
junc_idx = np.where((prof_diff > P_bg) & (depths_nm > 50))[0]
if len(junc_idx) > 0:
    xj = depths_nm[junc_idx[-1]]
    print(f"\n接合深さ xj ≈ {xj:.0f} nm  (B = P 背景 {P_bg:.0e} cm^-3)")

print(f"拡散後ピーク: {prof_diff.max():.3e} cm^-3  at {depths_nm[np.argmax(prof_diff)]:.1f} nm")

# --------------------------------------------------------------------------
# 5. VTU 保存 (ParaView 用)
# --------------------------------------------------------------------------
out_vtu = os.path.join(os.path.dirname(__file__), 'boron_sim.vtu')
cp.save_vtu(out_vtu, st)
print(f"\nVTK ファイル保存: {out_vtu}")

print("\n=== シミュレーション完了 ===")
