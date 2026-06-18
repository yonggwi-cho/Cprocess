#!/usr/bin/env python3
"""
Comprehensive tests for cprocess.Simulation covering all major features.

Run from the python/ directory (so the package resolves):

    cp ../build/_cprocess*.so cprocess/
    python3 test_comprehensive.py
"""

import os
import sys
import tempfile

import numpy as np

import cprocess as cp


def check(cond, msg):
    if not cond:
        raise AssertionError(msg)
    print(f"  ok: {msg}")


# ---------------------------------------------------------------------------
def test_bc_diffuse():
    """Set a Dirichlet BC, diffuse, verify the field stays finite."""
    print("test_bc_diffuse")
    sim = cp.Simulation()
    sim.mesh(x=0.4, y=0.4, z=0.8, nx=4, ny=4, nz=8)
    sim.region("silicon")
    sim.init("P", 1e15)
    sim.bc("P", "zmax", 1e20)
    sim.diffuse(time=5, temp=1000)
    P = sim.field("P")
    check(np.all(np.isfinite(P)), "P field finite after BC diffuse")
    check(P.max() > 0, "P > 0 after BC diffuse")


# ---------------------------------------------------------------------------
def test_dose():
    """Create a Gaussian implant and verify dose() > 0."""
    print("test_dose")
    sim = cp.Simulation()
    sim.mesh(x=0.4, y=0.4, z=0.8, nx=4, ny=4, nz=16)
    sim.region("silicon")
    sim.implant("B", dose=1e13, rp=0.05, drp=0.02)
    d = sim.dose("B")
    print(f"  dose = {d:.3e} cm^-2")
    check(d > 0, "dose() > 0")


# ---------------------------------------------------------------------------
def test_set_field():
    """set_field() round-trips back through field()."""
    print("test_set_field")
    sim = cp.Simulation()
    sim.mesh(x=0.4, y=0.4, z=0.8, nx=4, ny=4, nz=8)
    sim.region("silicon")
    n = sim.n_cells
    vals = np.linspace(1e14, 1e16, n)
    sim.set_field("Ga", vals)
    got = sim.field("Ga")
    check(np.allclose(got, vals), "set_field/field round-trip matches")


# ---------------------------------------------------------------------------
def test_analytic_window():
    """Implant with window=(0,0.2,0,0.4) concentrates dose inside window."""
    print("test_analytic_window")
    sim = cp.Simulation()
    sim.mesh(x=0.4, y=0.4, z=0.8, nx=6, ny=6, nz=16)
    sim.region("silicon")
    sim.implant("B", dose=1e13, rp=0.08, drp=0.03,
                window=(0, 0.2, 0, 0.4))
    xyz = sim.cell_centroids   # micrometres
    vol = sim.cell_volumes
    B = sim.field("B")
    inside = (xyz[:, 0] < 0.2) & (xyz[:, 1] < 0.4)
    q_in  = float(np.sum(B[inside]  * vol[inside]))
    q_out = float(np.sum(B[~inside] * vol[~inside]))
    print(f"  in={q_in:.3e}  out={q_out:.3e}  ratio={q_in/(q_out+1e-30):.1f}x")
    check(q_in > 0, "dose inside window > 0")
    check(q_in > q_out, "more dose inside window than outside")


# ---------------------------------------------------------------------------
def test_multi_species():
    """Init B and P, diffuse both, verify both fields are present and finite."""
    print("test_multi_species")
    sim = cp.Simulation()
    sim.mesh(x=0.4, y=0.4, z=0.8, nx=4, ny=4, nz=8)
    sim.region("silicon")
    sim.init("B", 1e15)
    sim.init("P", 5e14)
    sim.diffuse(time=5, temp=950)
    check("B" in sim.field_names(), "B field present")
    check("P" in sim.field_names(), "P field present")
    check(np.all(np.isfinite(sim.field("B"))), "B finite after diffuse")
    check(np.all(np.isfinite(sim.field("P"))), "P finite after diffuse")


# ---------------------------------------------------------------------------
def test_save_vtu():
    """save() writes a non-empty .vtu file."""
    print("test_save_vtu")
    sim = cp.Simulation()
    sim.mesh(x=0.4, y=0.4, z=0.8, nx=4, ny=4, nz=8)
    sim.region("silicon")
    sim.init("B", 1e15)
    with tempfile.NamedTemporaryFile(suffix=".vtu", delete=False) as f:
        path = f.name
    try:
        sim.save(path)
        check(os.path.exists(path), ".vtu file created")
        check(os.path.getsize(path) > 0, ".vtu file is non-empty")
    finally:
        if os.path.exists(path):
            os.unlink(path)


# ---------------------------------------------------------------------------
def test_clear_bc():
    """add bc, clear_bc, diffuse — must not error."""
    print("test_clear_bc")
    sim = cp.Simulation()
    sim.mesh(x=0.4, y=0.4, z=0.8, nx=4, ny=4, nz=8)
    sim.region("silicon")
    sim.init("B", 1e15)
    sim.implant("B", dose=1e13, rp=0.1, drp=0.03)  # non-trivial profile
    sim.bc("B", "zmax", 1e20)
    sim.clear_bc()
    sim.diffuse(time=5, temp=1000)
    B = sim.field("B")
    check(np.all(np.isfinite(B)), "B finite after clear_bc + diffuse")


# ---------------------------------------------------------------------------
def test_resist_diffuse():
    """photo → mask → implant → strip → diffuse → field finite."""
    print("test_resist_diffuse")
    sim = cp.Simulation()
    sim.mesh(x=1.0, y=1.0, z=0.5, nx=8, ny=8, nz=4)
    sim.region("silicon")
    sim.init("B", 1e15)
    sim.photo(resist=0.4)
    sim.mask(x1=0.0, x2=0.35)
    sim.mask(x1=0.65, x2=1.0)
    sim.implant("P", dose=5e14, energy=30, mc=True, ions=30000, seed=17)
    sim.strip()
    sim.diffuse(time=5, temp=1000)
    check("P" in sim.field_names(), "P field present after resist flow")
    P = sim.field("P")
    check(np.all(np.isfinite(P)), "P finite after resist diffuse")


# ---------------------------------------------------------------------------
def test_verbose_mode():
    """Simulation(verbose=True) accumulates a non-empty log."""
    print("test_verbose_mode")
    sim = cp.Simulation(verbose=True)
    sim.mesh(x=0.4, y=0.4, z=0.8, nx=4, ny=4, nz=8)
    sim.region("silicon")
    sim.init("B", 1e15)
    sim.implant("B", dose=1e13, rp=0.1, drp=0.03)  # non-trivial gradient
    sim.diffuse(time=5, temp=1000)
    check(len(sim.log) > 0, "verbose log is non-empty")


# ---------------------------------------------------------------------------
if __name__ == "__main__":
    test_bc_diffuse()
    test_dose()
    test_set_field()
    test_analytic_window()
    test_multi_species()
    test_save_vtu()
    test_clear_bc()
    test_resist_diffuse()
    test_verbose_mode()
    print("\nall comprehensive Simulation tests passed")
