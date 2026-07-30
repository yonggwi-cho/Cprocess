#!/usr/bin/env python3
"""
円形コンタクトホールのイオン注入サンプル。

1 µm × 1 µm のシリコン基板上に 0.4 µm の光レジストを形成し、
中心 (0.5, 0.5)、半径 0.2 µm の円形ホールを開口して
リン (P) を注入する。

mask_polygon() に 64 頂点の円形多角形を渡すことで円孔を近似する。

実行方法::

    cmake --build build --target _cprocess
    cp build/_cprocess*.so python/cprocess/
    PYTHONPATH=python python3 examples/circular_contact.py
"""

import math
import numpy as np
import cprocess as cp


def circle_polygon(cx: float, cy: float, r: float, n: int = 64):
    """中心 (cx, cy)、半径 r の円を n 頂点の多角形で近似する（µm 単位）。"""
    return [(cx + r * math.cos(2 * math.pi * i / n),
             cy + r * math.sin(2 * math.pi * i / n))
            for i in range(n)]


def main():
    sim = cp.Simulation(verbose=True)

    # 1 µm × 1 µm × 0.5 µm のシリコン基板
    sim.mesh(x=1.0, y=1.0, z=0.5, nx=10, ny=10, nz=5)
    sim.region("silicon")
    sim.init("B", 1e15)          # p 型バックグラウンド

    # 全面 0.4 µm の光レジストを形成
    sim.photo(resist=0.4)

    # 中心 (0.5, 0.5)、半径 0.2 µm の円形ホールを開口
    hole = circle_polygon(cx=0.5, cy=0.5, r=0.2, n=64)
    sim.mask_polygon(hole)

    # リンをイオン注入（レジスト下はブロック）
    result = sim.implant("P", dose=5e15, energy=30,
                         mc=True, ions=60000, seed=7)
    print(f"\nP implant: {result}")

    # レジストを除去
    sim.strip()

    # アニール
    sim.diffuse(time=20, temp=1000)

    # ドーズ量の空間分布を確認
    xyz = sim.cell_centroids      # shape (n_cells, 3)  [µm]
    vol = sim.cell_volumes        # [cm³]
    P   = sim.field("P")         # [cm⁻³]

    # 表面直下 (z < 0.15 µm) の面積積分ドーズ量で円孔内外を比較
    cx, cy, r = 0.5, 0.5, 0.2
    dist2 = (xyz[:, 0] - cx)**2 + (xyz[:, 1] - cy)**2
    near_surface = xyz[:, 2] < 0.15
    inside  = near_surface & (dist2 <  r**2)
    outside = near_surface & (dist2 >= r**2)

    area_hole = math.pi * (r * cp.UM)**2           # cm²
    bb = sim.bbox()  # µm (W-3 fix)
    area_total = ((bb[1][0]-bb[0][0]) * cp.UM) * ((bb[1][1]-bb[0][1]) * cp.UM)  # cm²
    area_mask  = area_total - area_hole

    dose_in  = float(np.sum(P[inside]  * vol[inside]))  / area_hole  if inside.any()  else 0
    dose_out = float(np.sum(P[outside] * vol[outside])) / area_mask  if outside.any() else 0

    print(f"\nリン 面積積分ドーズ量 (表面 0.15 µm 以内):")
    print(f"  円孔内 (open)  : {dose_in:.3e} cm⁻²")
    print(f"  レジスト下     : {dose_out:.3e} cm⁻²")
    ratio = dose_in / (dose_out + 1e-30)
    print(f"  ブロック比     : {ratio:.1f}× (open/masked)")

    sim.save("circular_contact.vtu")
    print("\ncircular_contact.vtu を書き出しました — ParaView で確認できます")


if __name__ == "__main__":
    main()
