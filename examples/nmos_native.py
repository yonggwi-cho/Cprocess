#!/usr/bin/env python3
"""
Native-Python NMOS source/drain flow.

Mirrors the C++ integration example but drives the simulator entirely through
the native `cprocess.Simulation` API — no text deck. Demonstrates a physical
photoresist gate: phosphorus is blocked over the channel and implanted into the
open source/drain windows.

Run from the repo root after building the Python module::

    cmake --build build --target _cprocess
    cp build/_cprocess*.so python/cprocess/
    PYTHONPATH=python python3 examples/nmos_native.py
"""

import numpy as np

import cprocess as cp


def main() -> None:
    sim = cp.Simulation(verbose=True)

    # 1 um x 1 um x 0.5 um silicon substrate.
    sim.mesh(x=1.0, y=1.0, z=0.5, nx=8, ny=8, nz=4)
    sim.region("silicon")
    sim.init("B", 1e15)  # p-type background

    # Channel-stop boron implant (blanket).
    sim.implant("B", dose=5e12, energy=40, mc=True, ions=40000,
                tilt=7, channeling=True, seed=42)

    # Photoresist gate: cover the channel, open source/drain.
    sim.photo(resist=0.4)
    sim.mask(x1=0.0, x2=0.35)    # source window
    sim.mask(x1=0.65, x2=1.0)    # drain window
    sd = sim.implant("P", dose=5e15, energy=30, mc=True, ions=80000, seed=7)
    sim.strip()
    print(f"\nS/D implant: {sd}")

    # Source/drain anneal.
    sim.diffuse(time=30, temp=1000)

    # Quantify masking: phosphorus dose under the gate vs in the openings.
    xyz = sim.cell_centroids                 # micrometres
    vol = sim.cell_volumes
    P = sim.field("P")
    gate = (xyz[:, 0] > 0.35) & (xyz[:, 0] < 0.65)
    q_gate = float(np.sum(P[gate] * vol[gate]))
    q_open = float(np.sum(P[~gate] * vol[~gate]))
    print(f"\nphosphorus in silicon:")
    print(f"  open S/D : {q_open:.3e} atoms")
    print(f"  gate     : {q_gate:.3e} atoms")
    print(f"  blocking : {q_open / (q_gate + 1e-30):.1f}x")

    sim.save("nmos_native.vtu")
    print("\nwrote nmos_native.vtu — open in ParaView")


if __name__ == "__main__":
    main()
