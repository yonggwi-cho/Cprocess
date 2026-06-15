"""
3D ブラウザ可視化スクリプト
==============================
シミュレーション結果を Plotly で自己完結型 HTML に変換する。
生成した HTML をブラウザで開くだけで 3D 可視化できる。

  python3 examples/visualize_3d.py
  → examples/boron_3d_viewer.html
"""

import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'python'))

import numpy as np
from scipy.interpolate import LinearNDInterpolator
import plotly.graph_objects as go
from plotly.subplots import make_subplots
import cprocess as cp

# --------------------------------------------------------------------------
# 1. シミュレーション実行 (前回と同じ条件)
# --------------------------------------------------------------------------
print("シミュレーション実行中...")
st = cp.SimState()
cp.run_deck("""
mesh box xmax=0.5um ymax=0.5um zmax=1.2um nx=10 ny=10 nz=80
init species=P conc=1e15
implant species=B energy=50keV dose=1e13 method=mc ions=80000 seed=7
diffuse time=30min temp=1000C dt=30s
""", st)

# --------------------------------------------------------------------------
# 2. メッシュデータ取得
# --------------------------------------------------------------------------
xyz = st.mesh.cell_cent          # (n,3) [cm]
B   = np.maximum(st.get_field('B'),  1e6)   # cm^-3 (ゼロ回避)
P   = np.maximum(st.get_field('P'),  1e6)

um = 1e-4   # cm → µm 変換

x_um = xyz[:, 0] * 1e4
y_um = xyz[:, 1] * 1e4
z_um = xyz[:, 2] * 1e4          # 底面=0, 表面=zmax
z_top_um = z_um.max()

depth_nm = (z_top_um - z_um) * 1e3   # 表面からの深さ [nm]

logB = np.log10(B)
logP = np.log10(P)

# --------------------------------------------------------------------------
# 3. 正規グリッドに補間 (等値面 / スライスに必要)
# --------------------------------------------------------------------------
print("正規グリッドに補間中...")
NX, NY, NZ = 30, 30, 80    # グリッド解像度
xi = np.linspace(x_um.min(), x_um.max(), NX)
yi = np.linspace(y_um.min(), y_um.max(), NY)
zi = np.linspace(z_um.min(), z_um.max(), NZ)
Xi, Yi, Zi = np.meshgrid(xi, yi, zi, indexing='ij')

pts  = np.column_stack([x_um, y_um, z_um])
interp_B = LinearNDInterpolator(pts, logB, fill_value=6.0)
interp_P = LinearNDInterpolator(pts, logP, fill_value=15.0)

BG = interp_B(Xi, Yi, Zi)  # (NX, NY, NZ)
PG = interp_P(Xi, Yi, Zi)

# 接合面 (B = P) の差分でゼロ交差を探す
JG = BG - PG                # >0 で B優勢, <0 で P優勢

# --------------------------------------------------------------------------
# 4. Plotly 図の構築
# --------------------------------------------------------------------------
print("HTML 生成中...")

fig = make_subplots(
    rows=1, cols=2,
    column_widths=[0.6, 0.4],
    specs=[[{"type": "scene"}, {"type": "xy"}]],
    subplot_titles=[
        "3D: B 濃度等値面 + 接合面 (B = P 背景)",
        "深さプロファイル (中心 XZ 断面)"
    ]
)

# ── 3D 等値面: ボロン (log B = 16, 17, 17.5) ──────────────────────────
opacity_levels = [(16.0, 0.08, "B=1e16"), (17.0, 0.15, "B=1e17"), (17.5, 0.25, "B=3e17")]
colors = ["#4fc3f7", "#0288d1", "#01579b"]
for (lev, opac, name), col in zip(opacity_levels, colors):
    fig.add_trace(go.Isosurface(
        x=Xi.ravel(), y=Yi.ravel(), z=Zi.ravel(),
        value=BG.ravel(),
        isomin=lev, isomax=lev,
        surface_count=1,
        colorscale=[[0, col], [1, col]],
        showscale=False,
        opacity=opac,
        caps=dict(x_show=False, y_show=False, z_show=False),
        name=name,
        showlegend=True,
    ), row=1, col=1)

# ── 3D 等値面: 接合面 B=P ──────────────────────────────────────────────
fig.add_trace(go.Isosurface(
    x=Xi.ravel(), y=Yi.ravel(), z=Zi.ravel(),
    value=JG.ravel(),
    isomin=0.0, isomax=0.0,
    surface_count=1,
    colorscale=[[0, "#ff7043"], [1, "#ff7043"]],
    showscale=False,
    opacity=0.5,
    caps=dict(x_show=False, y_show=False, z_show=False),
    name="接合面 (B=P)",
    showlegend=True,
), row=1, col=1)

# ── 2D 深さプロファイル (y=中心の XZ スライス) ───────────────────────
iy_mid = NY // 2
slice_z  = zi                        # µm
slice_logB = BG[:, iy_mid, :].mean(axis=0)   # x 平均
slice_logP = PG[:, iy_mid, :].mean(axis=0)
depth_slice_nm = (z_top_um - slice_z) * 1e3  # nm

# 深さプロファイル線グラフ
fig.add_trace(go.Scatter(
    x=depth_slice_nm, y=slice_logB,
    mode='lines', name='B (拡散後)',
    line=dict(color='#0288d1', width=2.5),
), row=1, col=2)
fig.add_trace(go.Scatter(
    x=depth_slice_nm, y=slice_logP,
    mode='lines', name='P 背景',
    line=dict(color='#e53935', width=2, dash='dash'),
), row=1, col=2)

# 接合深さの縦線 (shapes で追加)
junc_idx_arr = np.where(np.diff(np.sign(slice_logB - slice_logP)))[0]
for ji in junc_idx_arr:
    xj = float(depth_slice_nm[ji])
    fig.add_shape(
        type="line", xref="x2", yref="paper",
        x0=xj, x1=xj, y0=0, y1=1,
        line=dict(color="#ff7043", width=1.5, dash="dot"),
    )
    fig.add_annotation(
        xref="x2", yref="paper",
        x=xj, y=1.02, text=f"xj≈{xj:.0f}nm",
        showarrow=False, font=dict(size=11, color="#ff7043"),
    )

# --------------------------------------------------------------------------
# 5. レイアウト設定
# --------------------------------------------------------------------------
fig.update_layout(
    title=dict(
        text="<b>Cprocess: B 注入 (50 keV, 1e13 cm⁻²) + 熱拡散 (1000°C, 30 min) in Si</b>",
        font=dict(size=16),
    ),
    legend=dict(x=0.01, y=0.99, bgcolor="rgba(255,255,255,0.8)"),
    height=700,
    paper_bgcolor="white",
)

fig.update_scenes(
    xaxis_title="x [µm]",
    yaxis_title="y [µm]",
    zaxis_title="z [µm] (表面=" + f"{z_top_um:.2f}" + ")",
    camera=dict(eye=dict(x=1.6, y=1.6, z=0.8)),
    aspectmode="data",
    bgcolor="#f8f8f8",
)

fig.update_xaxes(title_text="表面からの深さ [nm]", row=1, col=2)
fig.update_yaxes(title_text="log₁₀(濃度 [cm⁻³])", row=1, col=2, range=[8, 19])
fig.add_shape(
    type="rect", xref="x2", yref="y2",
    x0=depth_slice_nm.min(), x1=depth_slice_nm.max(),
    y0=14.9, y1=15.1,
    fillcolor="#ffcdd2", opacity=0.3, line_width=0,
)

# --------------------------------------------------------------------------
# 6. HTML 保存
# --------------------------------------------------------------------------
out_path = os.path.join(os.path.dirname(__file__), 'boron_3d_viewer.html')
fig.write_html(out_path, include_plotlyjs='cdn', full_html=True)
print(f"\n✓ ブラウザで開いてください: {out_path}")
print("  (ファイルを直接ブラウザにドラッグ＆ドロップ)")
