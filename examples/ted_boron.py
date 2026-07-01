#!/usr/bin/env python3
"""
Transient Enhanced Diffusion (TED) of a shallow boron implant.

Ion-implant damage creates a large excess of silicon self-interstitials. During
the initial seconds-to-minutes of the anneal these interstitials supersaturate
the lattice and dramatically accelerate boron diffusion (boron moves by pairing
with interstitials). The enhancement decays as the interstitials diffuse to the
surface sink and recombine — so the profile broadens fast at first, then
settles into ordinary equilibrium diffusion.

This script compares a shallow B profile annealed with and without TED and
prints the resulting junction broadening.

Run from the repo root after building::

    cmake --build build --target _cprocess
    cp build/_cprocess*.so python/cprocess/
    PYTHONPATH=python python3 examples/ted_boron.py
"""

import numpy as np
import cprocess as cp


def spread_um(sim, species):
    """Mass-weighted standard deviation of the profile along depth (um)."""
    c = sim.field(species)
    v = sim.cell_volumes
    z = sim.cell_centroids[:, 2]
    m = float(np.sum(c * v))
    mean = float(np.sum(c * v * z)) / m
    var = float(np.sum(c * v * z * z)) / m - mean * mean
    return float(np.sqrt(max(0.0, var)))


def anneal(ted: bool):
    sim = cp.Simulation()
    sim.mesh(x=0.3, y=0.3, z=1.0, nx=6, ny=6, nz=60)
    sim.region("silicon")
    # Shallow 10 keV-class boron implant, Rp = 50 nm.
    sim.implant("B", dose=1e14, rp=0.05, drp=0.02, damage=ted)
    s0 = spread_um(sim, "B")
    # 60 s spike anneal at 950 C.
    sim.diffuse(time=1.0, temp=950, ted=ted)
    s1 = spread_um(sim, "B")
    return sim, s0, s1


def main():
    sim_eq,  s0, s_eq  = anneal(ted=False)
    sim_ted, _,  s_ted = anneal(ted=True)

    print("shallow boron, 60 s @ 950 C")
    print(f"  as-implanted spread : {s0 * 1000:.1f} nm")
    print(f"  equilibrium anneal  : {s_eq * 1000:.1f} nm  (+{(s_eq-s0)*1000:.1f} nm)")
    print(f"  TED anneal          : {s_ted * 1000:.1f} nm  (+{(s_ted-s0)*1000:.1f} nm)")
    print(f"  TED / equilibrium broadening : {(s_ted-s0)/(s_eq-s0+1e-30):.1f}x")

    Ipeak = sim_ted.field("I").max()
    print(f"  seeded interstitial peak (residual after anneal): {Ipeak:.2e} cm^-3")

    sim_ted.save("ted_boron.vtu")
    print("\nwrote ted_boron.vtu (fields: B, I) — open in ParaView")


if __name__ == "__main__":
    main()
