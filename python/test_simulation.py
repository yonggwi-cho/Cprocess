#!/usr/bin/env python3
"""
Smoke + regression tests for the native cprocess.Simulation API.

Run from the python/ directory (so the package resolves):

    cp ../build/_cprocess*.so cprocess/
    python3 test_simulation.py
"""

import numpy as np

import cprocess as cp


def check(cond, msg):
    if not cond:
        raise AssertionError(msg)
    print(f"  ok: {msg}")


def test_basic_flow():
    print("test_basic_flow")
    sim = cp.Simulation()
    sim.mesh(x=0.4, y=0.4, z=0.8, nx=6, ny=6, nz=12)
    sim.region("silicon")
    sim.init("P", 1e15)
    check(sim.n_cells > 0, "mesh created")
    check("P" in sim.field_names(), "P field present after init")

    r = sim.implant("B", dose=1e13, energy=50, mc=True, ions=20000, seed=3)
    check(r.deposited > 0, "MC implant deposited ions")
    check(r.rp_nm > 0, "MC implant reports Rp")

    sim.diffuse(time=10, temp=1000)
    B = sim.field("B")
    check(np.all(np.isfinite(B)), "B field finite after diffuse")
    check(B.max() > 0, "B present after diffuse")


def test_physical_resist():
    print("test_physical_resist")
    sim = cp.Simulation()
    sim.mesh(x=1.0, y=1.0, z=0.5, nx=8, ny=8, nz=4)
    sim.region("silicon")
    sim.photo(resist=0.4)
    sim.mask(x1=0.0, x2=0.35)
    sim.mask(x1=0.65, x2=1.0)
    sim.implant("P", dose=5e15, energy=30, mc=True, ions=40000, seed=17)
    sim.strip()

    xyz = sim.cell_centroids
    vol = sim.cell_volumes
    P = sim.field("P")
    gate = (xyz[:, 0] > 0.35) & (xyz[:, 0] < 0.65)
    q_gate = float(np.sum(P[gate] * vol[gate]))
    q_open = float(np.sum(P[~gate] * vol[~gate]))
    print(f"  open={q_open:.3e} gate={q_gate:.3e} ratio={q_open/(q_gate+1e-30):.1f}x")
    check(q_open > 0, "open S/D received phosphorus")
    check(q_open > 3.0 * q_gate, "resist blocks ions under the gate")


def test_gaussian_implant():
    print("test_gaussian_implant")
    sim = cp.Simulation()
    sim.mesh(x=0.4, y=0.4, z=0.8, nx=4, ny=4, nz=20)
    sim.region("silicon")
    atoms = sim.implant("B", dose=1e13, rp=0.1, drp=0.03)  # analytic Gaussian
    check(atoms > 0, "Gaussian implant returns atom count")
    check(sim.field("B").max() > 0, "B profile present")


def test_chaining():
    print("test_chaining")
    sim = cp.Simulation()
    (sim.mesh(x=0.4, y=0.4, z=0.8, nx=4, ny=4, nz=8)
        .region("silicon")
        .init("As", 1e15)
        .diffuse(time=5, temp=950))
    check("As" in sim.field_names(), "method chaining works")


if __name__ == "__main__":
    test_basic_flow()
    test_physical_resist()
    test_gaussian_implant()
    test_chaining()
    print("\nall native Simulation tests passed")
