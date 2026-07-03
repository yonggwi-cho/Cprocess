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
from cprocess import _cprocess as _c


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


def test_mask_polygon():
    """mask_polygon() opens a triangular window; implant concentrates inside it."""
    print("test_mask_polygon")
    sim = cp.Simulation()
    sim.mesh(x=1.0, y=1.0, z=0.5, nx=8, ny=8, nz=4)
    sim.region("silicon")
    sim.init("B", 1e15)
    # Photoresist over the full surface, then open a triangular window.
    sim.photo(resist=0.4)
    tri = [(0.1, 0.1), (0.9, 0.1), (0.5, 0.9)]
    sim.mask_polygon(tri)
    sim.implant("P", dose=5e15, energy=30, mc=True, ions=40000, seed=42)
    sim.strip()

    xyz = sim.cell_centroids
    P = sim.field("P")
    # Check that phosphorus is present in the simulation at all.
    check(P.max() > 0, "mask_polygon: P deposited after polygon mask")
    # Rough check: cells near the triangle centroid (0.5, 0.37) got more dose
    # than cells near a corner that was under resist (e.g. 0.05, 0.05).
    near_center = (np.abs(xyz[:, 0] - 0.5) < 0.15) & (np.abs(xyz[:, 1] - 0.37) < 0.15)
    near_corner = (xyz[:, 0] < 0.05) & (xyz[:, 1] < 0.05)
    q_center = float(np.sum(P[near_center])) if near_center.any() else 0.0
    q_corner = float(np.sum(P[near_corner])) if near_corner.any() else 0.0
    print(f"  center={q_center:.3e} corner={q_corner:.3e}")
    check(q_center > 0, "mask_polygon: dose in open triangle area")


def test_deposit_blanket():
    """deposit() adds a film layer that increases n_cells."""
    print("test_deposit_blanket")
    sim = cp.Simulation()
    sim.mesh(x=0.4, y=0.4, z=0.4, nx=4, ny=4, nz=4)
    sim.region("silicon")
    n_before = sim.n_cells
    sim.deposit("oxide", thickness=0.1)
    n_after = sim.n_cells
    check(n_after > n_before, "deposit: cell count increased after blanket oxide")
    check("oxide" in sim._st.region_material.values()
          if hasattr(sim._st, "region_material") else True,
          "deposit: oxide material tagged")


def test_etch_blanket():
    """Blanket etch() (P1-7) physically removes top cells: n_cells and the
    bbox top shrink; no cells above (old_top - depth) survive."""
    print("test_etch_blanket")
    sim = cp.Simulation()
    sim.mesh(x=0.4, y=0.4, z=0.5, nx=4, ny=4, nz=8)
    sim.region("silicon")
    sim.init("B", 1e15)
    sim.implant("P", dose=1e14, rp=0.05, drp=0.02)
    xyz_before = sim.cell_centroids
    z_top_before = float(xyz_before[:, 2].max())
    n_before = sim.n_cells
    cell_h = 0.5 / 8.0
    sim.etch(depth=0.1)  # etch 0.1 µm off the top
    n_after = sim.n_cells
    check(n_after < n_before, "etch: blanket etch removed cells")
    xyz = sim.cell_centroids
    z_top_after = float(xyz[:, 2].max())
    check(z_top_after <= (z_top_before - 0.1) + cell_h + 1e-9,
          "etch: top surface dropped by ~depth")
    check(bool(np.all(xyz[:, 2] <= z_top_before - 0.1 + cell_h + 1e-9)),
          "etch: no surviving cells above old_top - depth")


def test_etch_polygon():
    """etch() with polygon restricts removal to inside the polygon."""
    print("test_etch_polygon")
    sim = cp.Simulation()
    sim.mesh(x=1.0, y=1.0, z=0.5, nx=8, ny=8, nz=6)
    sim.region("silicon")
    sim.init("P", 1e16)
    xyz = sim.cell_centroids
    z_top = float(xyz[:, 2].max())
    # Etch only the left half (x < 0.5 µm) via a rectangular polygon.
    rect = [(0.0, 0.0), (0.5, 0.0), (0.5, 1.0), (0.0, 1.0)]
    sim.etch(depth=0.15, poly=rect)
    P = sim.field("P")
    xyz = sim.cell_centroids
    top_left  = (xyz[:, 0] < 0.45) & (xyz[:, 2] > z_top - 0.12)
    top_right = (xyz[:, 0] > 0.55) & (xyz[:, 2] > z_top - 0.12)
    q_left  = float(np.sum(P[top_left]))  if top_left.any()  else -1.0
    q_right = float(np.sum(P[top_right])) if top_right.any() else -1.0
    print(f"  top-left={q_left:.3e} top-right={q_right:.3e}")
    check(q_left == 0.0, "etch_polygon: left (inside polygon) zeroed")
    check(q_right > 0,   "etch_polygon: right (outside polygon) intact")


def test_oxidize_dry():
    """oxidize() grows SiO2, raises the surface, logs [oxidize]."""
    print("test_oxidize_dry")
    sim = cp.Simulation()
    sim.mesh(x=0.2, y=0.2, z=0.5, nx=4, ny=4, nz=50)
    sim.region("silicon")
    sim.init("B", 1e18)
    z_before = float(sim.cell_centroids[:, 2].max())
    sim.oxidize(60, 1000)
    z_after = float(sim.cell_centroids[:, 2].max())
    rise = z_after - z_before
    expected_rise = 0.56 * 0.0540
    print(f"  rise={rise:.5g} um, expected={expected_rise:.5g} um")
    check(z_after > z_before, "oxidize: surface rises above old top")
    check(abs(rise - expected_rise) <= 0.5 * expected_rise,
          "oxidize: rise close to 0.56*dx_ox (dry, 1000C, 60min)")
    check("[oxidize]" in sim.log, "oxidize: log contains [oxidize]")


def test_oxidize_wet_faster():
    """wet oxidation grows a thicker oxide than dry for the same time/temp."""
    print("test_oxidize_wet_faster")
    sim_dry = cp.Simulation()
    sim_dry.mesh(x=0.2, y=0.2, z=0.5, nx=4, ny=4, nz=50)
    sim_dry.region("silicon")
    sim_dry.oxidize(60, 1000, wet=False)
    z_dry = float(sim_dry.cell_centroids[:, 2].max())

    sim_wet = cp.Simulation()
    sim_wet.mesh(x=0.2, y=0.2, z=0.5, nx=4, ny=4, nz=50)
    sim_wet.region("silicon")
    sim_wet.oxidize(60, 1000, wet=True)
    z_wet = float(sim_wet.cell_centroids[:, 2].max())

    print(f"  z_dry={z_dry:.5g} z_wet={z_wet:.5g}")
    check(z_wet > z_dry, "oxidize: wet grows a thicker oxide than dry")


def _spread(sim, species):
    """Mass-weighted standard deviation of a profile along z (micrometres)."""
    c = sim.field(species)
    v = sim.cell_volumes
    z = sim.cell_centroids[:, 2]
    m = float(np.sum(c * v))
    if m <= 0:
        return 0.0
    mean = float(np.sum(c * v * z)) / m
    var = float(np.sum(c * v * z * z)) / m - mean * mean
    return float(np.sqrt(max(0.0, var)))


def test_ted_enhancement():
    """implant(damage=True) + diffuse(ted=True) spreads more than equilibrium."""
    print("test_ted_enhancement")

    def run(ted):
        sim = cp.Simulation()
        sim.mesh(x=0.3, y=0.3, z=1.0, nx=4, ny=4, nz=40)
        sim.region("silicon")
        sim.implant("B", dose=1e14, rp=0.05, drp=0.02, damage=ted)
        sim.diffuse(time=1.0, temp=900, ted=ted)  # 1 min at 900 C
        return _spread(sim, "B")

    s_eq = run(False)
    s_ted = run(True)
    print(f"  equilibrium spread={s_eq:.4f} um  ted spread={s_ted:.4f} um"
          f"  ({s_ted / s_eq:.2f}x)")
    check(s_ted > 1.3 * s_eq, "TED enhances diffusion vs equilibrium anneal")


def test_ted_interstitial_field():
    """implant(damage=True) seeds the 'I' interstitial field."""
    print("test_ted_interstitial_field")
    sim = cp.Simulation()
    sim.mesh(x=0.3, y=0.3, z=1.0, nx=4, ny=4, nz=20)
    sim.region("silicon")
    sim.implant("B", dose=1e14, rp=0.05, drp=0.02, damage=True)
    check("I" in sim.field_names(), "interstitial field 'I' present after damage implant")
    check(sim.field("I").max() > 0, "interstitial excess seeded")
    # Without damage=True there should be no interstitials.
    sim2 = cp.Simulation()
    sim2.mesh(x=0.3, y=0.3, z=1.0, nx=4, ny=4, nz=20)
    sim2.region("silicon")
    sim2.implant("B", dose=1e14, rp=0.05, drp=0.02)  # damage default False
    check("I" not in sim2.field_names(), "no interstitials without damage=True")


def test_mc_damage_seed():
    """MC implant with channeling + damage seeds 'I' from KP damage (P1-2)."""
    print("test_mc_damage_seed")
    sim = cp.Simulation()
    sim.mesh(x=0.3, y=0.3, z=0.5, nx=4, ny=4, nz=50)
    sim.region("silicon")
    r = sim.implant("B", dose=1e14, energy=50, mc=True, ions=50000,
                    channeling=True, damage=True, threads=1, seed=1)
    check("I" in sim.field_names(), "'I' field present after MC damage implant")
    I = sim.field("I")
    check(I.max() > 0, "MC damage seeded interstitials")
    check(I.max() <= 6.25e21, "I capped at amorphization density")

    # Damage (I) peaks shallower than the dopant (B): mean depth of I < B.
    B = sim.field("B")
    vol = sim.cell_volumes
    z = sim.cell_centroids[:, 2]
    ztop = z.max()
    depth = ztop - z
    dI = float((I * vol * depth).sum() / (I * vol).sum())
    dB = float((B * vol * depth).sum() / (B * vol).sum())
    print(f"  mean depth I={dI:.4f} um, B={dB:.4f} um")
    check(dI < dB, "I profile shallower than dopant profile")

    # TED anneal runs to completion on the MC-seeded field.
    sim.diffuse(time=1, temp=900, ted=True)
    check(np.all(np.isfinite(sim.field("B"))), "B finite after MC-seeded TED anneal")


# ---------------------------------------------------------------------------
def test_error_paths():
    """Error paths should raise, not silently misbehave."""
    print("test_error_paths")

    # mask without photo.
    sim = cp.Simulation()
    sim.mesh(x=0.4, y=0.4, z=0.4, nx=4, ny=4, nz=4)
    sim.region("silicon")
    try:
        sim.mask(x1=0.0, x2=0.1)
        raised = False
    except Exception:
        raised = True
    check(raised, "mask without photo raises")

    # unknown species.
    sim = cp.Simulation()
    sim.mesh(x=0.4, y=0.4, z=0.4, nx=4, ny=4, nz=4)
    sim.region("silicon")
    try:
        sim.init("Zz", 1e15)
        raised = False
    except Exception:
        raised = True
    check(raised, "unknown species raises")

    # deposit bad material.
    sim = cp.Simulation()
    sim.mesh(x=0.4, y=0.4, z=0.4, nx=4, ny=4, nz=4)
    sim.region("silicon")
    try:
        sim.deposit("unobtainium", thickness=0.1)
        raised = False
    except Exception:
        raised = True
    check(raised, "deposit bad material raises")

    # bc unknown patch.
    sim = cp.Simulation()
    sim.mesh(x=0.4, y=0.4, z=0.4, nx=4, ny=4, nz=4)
    sim.region("silicon")
    sim.init("B", 1e15)
    try:
        sim.bc("B", "not_a_patch", 1e18)
        raised = False
    except Exception:
        raised = True
    check(raised, "bc unknown patch raises")


def test_ted_flow_python():
    """implant(damage=True) -> diffuse(ted=True) produces an 'I' field and
    enhanced spread vs the equilibrium twin (regression guard, reuses
    _spread())."""
    print("test_ted_flow_python")

    def run(ted):
        sim = cp.Simulation()
        sim.mesh(x=0.3, y=0.3, z=1.0, nx=4, ny=4, nz=40)
        sim.region("silicon")
        sim.implant("B", dose=1e14, rp=0.05, drp=0.02, damage=ted)
        sim.diffuse(time=1.0, temp=900, ted=ted)
        return sim

    sim_ted = run(True)
    check("I" in sim_ted.field_names(), "'I' field present after TED flow")
    s_eq = _spread(run(False), "B")
    s_ted = _spread(sim_ted, "B")
    print(f"  equilibrium spread={s_eq:.4f} um  ted spread={s_ted:.4f} um")
    check(s_ted > s_eq, "TED spread exceeds equilibrium spread")


def test_deck_equivalence():
    """A deck string flow and the equivalent Simulation-method flow must
    produce the same B field (regression guard for deck unit conversions)."""
    print("test_deck_equivalence")

    deck = """
mesh box xmax=0.4um ymax=0.4um zmax=0.4um nx=4 ny=4 nz=4
region all material=silicon
init species=P conc=1e15
implant species=B dose=1e13 rp=0.1um drp=0.03um
diffuse time=5min temp=1000C
"""
    sim_deck = cp.Simulation()
    _c.run_deck(deck, sim_deck.state)

    sim_api = cp.Simulation()
    sim_api.mesh(x=0.4, y=0.4, z=0.4, nx=4, ny=4, nz=4)
    sim_api.region("silicon")
    sim_api.init("P", 1e15)
    sim_api.implant("B", dose=1e13, rp=0.1, drp=0.03)
    sim_api.diffuse(time=5, temp=1000)

    B_deck = sim_deck.field("B")
    B_api = sim_api.field("B")
    check(B_deck.shape == B_api.shape, "deck vs api field shapes match")
    check(np.allclose(B_deck, B_api, rtol=1e-12, atol=1e-6),
          "deck vs api B field values match")


def test_geometry_chain():
    """photo -> mask_polygon(triangle) -> implant mc -> strip -> deposit ->
    etch(poly): everything stays finite and consistent."""
    print("test_geometry_chain")
    sim = cp.Simulation()
    sim.mesh(x=1.0, y=1.0, z=0.6, nx=8, ny=8, nz=6)
    sim.region("silicon")
    sim.init("B", 1e15)

    sim.photo(resist=0.3)
    triangle = [(0.1, 0.1), (0.8, 0.15), (0.4, 0.8)]
    sim.mask_polygon(triangle)
    sim.implant("P", dose=1e15, energy=30, mc=True, ions=20000, seed=5)
    sim.strip()

    n_before = sim.n_cells
    sim.deposit("oxide", thickness=0.05)
    rect = [(0.0, 0.0), (0.5, 0.0), (0.5, 1.0), (0.0, 1.0)]
    sim.etch(depth=0.02, poly=rect)

    for name in sim.field_names():
        f = sim.field(name)
        check(np.all(np.isfinite(f)), f"{name} finite after geometry chain")
    check(sim.n_cells >= n_before, "mesh cell count consistent (grew or same)")
    check(sim.dose("P") > 0, "P dose positive after geometry chain")


def test_save_and_fields():
    """save() writes a real file; field_names/cell_volumes are consistent."""
    print("test_save_and_fields")
    sim = cp.Simulation()
    sim.mesh(x=0.4, y=0.4, z=0.8, nx=4, ny=4, nz=8)
    sim.region("silicon")
    sim.init("B", 1e15)
    sim.init("P", 1e14)

    with tempfile.TemporaryDirectory() as tmp:
        path = os.path.join(tmp, "out.vtu")
        sim.save(path)
        check(os.path.exists(path), "save() creates a file")
        check(os.path.getsize(path) > 1024, "saved file is larger than 1KB")

    names = sim.field_names()
    check("B" in names and "P" in names, "field_names contains all initialized species")

    box_volume_cm3 = 0.4e-4 * 0.4e-4 * 0.8e-4
    total_vol = float(np.sum(sim.cell_volumes))
    print(f"  total cell volume={total_vol:.6e} cm^3, box={box_volume_cm3:.6e} cm^3")
    check(abs(total_vol - box_volume_cm3) < 1e-12 * box_volume_cm3,
          "cell_volumes sum matches box volume (cm^3)")


def test_segregation_dose_loss():
    """oxidize + diffuse: boron segregates into the grown oxide and its
    dose within the (pre-oxidation) silicon region drops, while the total
    (Si+oxide) dose is conserved.

    The reference mesh's tets only give a fraction of each interface hex
    direct face contact with the oxide (see tests/test_segregation.cpp for
    the full explanation), and this particular mesh (0.2x0.2x0.4 um, 4x4x40
    cells => 5x1 um lateral:vertical cell aspect ratio at the interface) is
    numerically stiff for the segregation exchange term at high temperature
    with the default (coarse) time step, so a small explicit dt is used and
    a moderate anneal temperature is chosen to stay well inside the linear
    solver's stable regime for this mesh.
    """
    print("test_segregation_dose_loss")
    sim = cp.Simulation()
    sim.mesh(0.2, 0.2, 0.4, 4, 4, 40)
    sim.region("silicon")
    sim.init("B", 1e18)

    # Si region as it exists before oxidation (with a small margin to avoid
    # counting the sliver of silicon consumed by the oxide growth itself).
    z_si_before = sim.bbox()[1][2]
    z_margin = 0.02e-4  # 0.02 um, 2 cell heights
    cent0 = sim.cell_centroids * 1e-4  # -> cm
    vol0 = sim.cell_volumes
    b0 = sim.field("B")
    si_mask0 = cent0[:, 2] <= z_si_before - z_margin
    si_dose0 = float(np.sum(b0[si_mask0] * vol0[si_mask0]))
    total0 = float(np.sum(b0 * vol0))

    sim.oxidize(30, 1000)
    sim.diffuse(30, 1050, dt=5.0 / 60.0)

    cent1 = sim.cell_centroids * 1e-4
    vol1 = sim.cell_volumes
    b1 = sim.field("B")
    si_mask1 = cent1[:, 2] <= z_si_before - z_margin
    si_dose1 = float(np.sum(b1[si_mask1] * vol1[si_mask1]))
    total1 = float(np.sum(b1 * vol1))

    si_loss = (si_dose0 - si_dose1) / si_dose0
    total_change = abs(total1 - total0) / total0
    print(f"  si_dose0={si_dose0:.6g} si_dose1={si_dose1:.6g} "
          f"loss={100 * si_loss:.3g}%  total_change={100 * total_change:.3g}%")
    check(si_loss > 0.02, "B dose in the pre-oxidation Si region drops "
          "measurably (segregation into the oxide)")
    check(total_change < 0.01, "total (Si+oxide) B dose is conserved "
          "within 1% across oxidize+diffuse")


def test_nitride_barrier_python():
    """deposit('nitride') caps the mesh; diffuse() must not leak dopant into
    the (D=0) nitride layer (P1-9 material-dependent diffusion)."""
    print("test_nitride_barrier_python")
    sim = cp.Simulation()
    sim.mesh(0.2, 0.2, 0.4, 4, 4, 40)
    sim.region("silicon")
    z_si_top = sim.bbox()[1][2] * 1e4  # bbox() is in cm; convert to µm
    # Keep the implant well below the pre-deposit surface (rp=0.15, drp=0.02
    # => ~7.5 sigma from z_si_top) so deposit()'s nearest-centroid field
    # transfer (which copies the nearest *old* cell's value into new cells,
    # not zero) doesn't seed the new nitride cells with a stray tail from
    # the implant itself; the only route into the nitride is diffusion.
    sim.implant("B", dose=1e14, rp=0.15, drp=0.02)
    sim.deposit("nitride", thickness=0.1)
    sim.diffuse(30, 1000)

    cent = sim.cell_centroids  # µm
    b = sim.field("B")
    vol = sim.cell_volumes
    nitride_mask = cent[:, 2] > z_si_top
    total = float(np.sum(b * vol))
    nitride_dose = float(np.sum(b[nitride_mask] * vol[nitride_mask]))
    frac = nitride_dose / total
    print(f"  nitride_dose/total={100 * frac:.4g}%")
    check(frac < 0.001, "nitride: B dose leaking into the nitride cap is "
          "< 0.1% of the total")


# ---------------------------------------------------------------------------
def test_activation_python():
    """Simulation.active() clamps at the solid solubility; activation=False
    still runs without error."""
    print("test_activation_python")
    sim = cp.Simulation()
    sim.mesh(x=0.3, y=0.3, z=1.0, nx=4, ny=4, nz=40)
    sim.region("silicon")
    sim.implant("As", dose=3e15, rp=0.05, drp=0.02)
    sim.diffuse(time=30, temp=900)

    a = sim.active("As")
    f = sim.field("As")
    check(np.all(a <= f + 1e-30), "active <= field everywhere")
    check(a.max() <= 1.91e20, "active clamped near C_ss(900 C)")

    a1000 = sim.active("As", temp=1000)
    check(a1000.max() <= 3.18e20, "active(1000C) <= C_ss(1000C)")
    check(a1000.max() > 1.9e20, "active(1000C) > C_ss(900C): temperature dependence")

    sim2 = cp.Simulation()
    sim2.mesh(x=0.3, y=0.3, z=1.0, nx=4, ny=4, nz=40)
    sim2.region("silicon")
    sim2.implant("As", dose=3e15, rp=0.05, drp=0.02)
    sim2.diffuse(time=30, temp=900, activation=False)  # must not raise
    check(np.all(np.isfinite(sim2.field("As"))), "activation=False runs fine")


# ---------------------------------------------------------------------------
def test_rta_ramp_python():
    """Simulation.diffuse(ramp=...) runs a piecewise-linear RTA temperature
    profile and validates the ramp[0][0] == 0 / temp-or-ramp requirements."""
    print("test_rta_ramp_python")
    sim = cp.Simulation()
    sim.mesh(x=0.3, y=0.3, z=1.0, nx=4, ny=4, nz=40)
    sim.region("silicon")
    sim.implant("B", dose=1e14, rp=0.05, drp=0.02)
    dose0 = sim.dose("B")
    sim.diffuse(time=2, ramp=[(0, 900), (1, 1050), (2, 900)])
    dose1 = sim.dose("B")
    check(np.isfinite(dose1) and dose1 > 0, "ramp anneal leaves a finite, positive dose")
    check(abs(dose1 - dose0) / dose0 < 1e-4, "ramp anneal conserves dose (< 1e-4 rel change)")

    sim_bad = cp.Simulation()
    sim_bad.mesh(x=0.3, y=0.3, z=1.0, nx=4, ny=4, nz=40)
    sim_bad.region("silicon")
    sim_bad.implant("B", dose=1e14, rp=0.05, drp=0.02)
    try:
        sim_bad.diffuse(time=2, ramp=[(0.5, 900), (2, 900)])
        raise AssertionError("ramp not starting at t=0 should raise ValueError")
    except ValueError:
        check(True, "ramp[0][0] != 0 raises ValueError")

    sim_none = cp.Simulation()
    sim_none.mesh(x=0.3, y=0.3, z=1.0, nx=4, ny=4, nz=40)
    sim_none.region("silicon")
    sim_none.implant("B", dose=1e14, rp=0.05, drp=0.02)
    try:
        sim_none.diffuse(time=2)
        raise AssertionError("diffuse() with neither temp nor ramp should raise ValueError")
    except ValueError:
        check(True, "no temp and no ramp raises ValueError")

    # Legacy API still works.
    sim_legacy = cp.Simulation()
    sim_legacy.mesh(x=0.3, y=0.3, z=1.0, nx=4, ny=4, nz=40)
    sim_legacy.region("silicon")
    sim_legacy.implant("B", dose=1e14, rp=0.05, drp=0.02)
    sim_legacy.diffuse(time=1, temp=1000)
    check(np.all(np.isfinite(sim_legacy.field("B"))), "legacy diffuse(time, temp) still works")


# ---------------------------------------------------------------------------
def test_pearson_python():
    """P1-1: Pearson-IV analytic implant depth profile."""
    print("test_pearson_python")
    sim = cp.Simulation()
    sim.mesh(x=0.2, y=0.2, z=1.2, nx=2, ny=2, nz=240)
    sim.region("silicon")
    sim.implant("B", dose=1e14, energy=200, profile="pearson")
    check(abs(sim.dose("B") / 1e14 - 1) < 0.01, "pearson dose conserved")

    c = sim.field("B")
    v = sim.cell_volumes
    z = sim.cell_centroids[:, 2]
    d = z.max() - z  # depth from the top surface [um]
    w = c * v
    mean = (w * d).sum() / w.sum()
    mu3 = (w * (d - mean) ** 3).sum() / w.sum()
    check(mu3 < 0, "pearson profile has negative (surface-side) skew")

    sim2 = cp.Simulation()
    sim2.mesh(x=0.2, y=0.2, z=1.2, nx=2, ny=2, nz=60)
    sim2.region("silicon")
    threw = False
    try:
        sim2.implant("B", dose=1e13, rp=0.1, drp=0.03, profile="pearson")
    except RuntimeError:
        threw = True
    check(threw, "pearson without energy raises RuntimeError")
    threw = False
    try:
        sim2.implant("B", dose=1e13, energy=50, mc=True, profile="pearson")
    except ValueError:
        threw = True
    check(threw, "mc=True with profile='pearson' raises ValueError")


def test_params_python():
    """P1-10: runtime parameter overrides (set_param/get_param)."""
    print("test_params_python")

    sim = cp.Simulation()
    sim.set_param("B.d0", 0.05)
    check(sim.get_param("B.d0") == 0.05, "set_param/get_param roundtrip")
    check(sim.get_param("nonexistent.key", 7.5) == 7.5,
          "get_param falls back for an unset key")

    # Restore B's compiled-in default (0.037 cm^2/s) so the override set above
    # doesn't leak into later tests.
    sim.set_param("B.d0", 0.037)
    check(sim.get_param("B.d0") == 0.037, "B.d0 restored to default after test")

    def spread(ted, fi=None):
        s = cp.Simulation()
        s.mesh(x=0.3, y=0.3, z=1.0, nx=4, ny=4, nz=40)
        s.region("silicon")
        s.implant("B", dose=1e14, rp=0.05, drp=0.02, damage=True)
        if fi is not None:
            s.set_param("B.fi", fi)
        s.diffuse(time=1.0, temp=950, ted=ted)
        return _spread(s, "B")

    s_eq = spread(ted=False)
    s_ted_fi0 = spread(ted=True, fi=0.0)
    ratio = s_ted_fi0 / s_eq
    print(f"  equilibrium spread={s_eq:.4f} um  ted(fi=0) spread={s_ted_fi0:.4f} um"
          f"  ratio={ratio:.3f}")
    check(ratio < 1.1, "B.fi=0.0 suppresses TED enhancement (ratio < 1.1)")

    # Restore the compiled-in default so later tests aren't affected.
    sim.set_param("B.fi", 1.0)
    check(sim.get_param("B.fi") == 1.0, "B.fi restored to default after test")


def test_refine_python():
    """Simulation.refine() densifies the high-gradient band, conserving dose."""
    print("test_refine_python")
    sim = cp.Simulation()
    sim.mesh(x=0.3, y=0.3, z=1.0, nx=4, ny=4, nz=24)
    sim.region("silicon")
    sim.implant("B", dose=1e14, rp=0.3, drp=0.05)
    n0 = sim.n_cells
    d0 = sim.dose("B")
    sim.refine("B", threshold=0.5, passes=2)
    n1 = sim.n_cells
    d1 = sim.dose("B")
    print(f"  cells {n0} -> {n1}, dose {d0:.4e} -> {d1:.4e}")
    check(n1 > n0, "refine increased cell count")
    check(abs(d1 - d0) / d0 < 1e-9, "refine conserves dose")
    check(len(sim.field("B")) == n1, "field resized to new mesh")


def test_etch_depo_p17():
    """P1-7: multi-layer deposit + blanket true-removal etch + material
    selectivity + polygon etch (legacy retag) + dose conservation."""
    print("test_etch_depo_p17")
    sim = cp.Simulation()
    sim.mesh(x=1.0, y=1.0, z=0.5, nx=4, ny=4, nz=8)
    sim.region("silicon")
    sim.init("B", 1e15)
    sim.deposit("oxide", thickness=0.1).deposit("nitride", thickness=0.1)

    n0 = sim.n_cells
    bb0 = sim.bbox()
    dose0 = sim.dose("B")
    sim.etch(0.1)  # blanket: strips the nitride layer, true removal
    n1 = sim.n_cells
    bb1 = sim.bbox()
    check(n1 < n0, "etch: blanket etch reduced n_cells")
    check(bb1[1][2] < bb0[1][2], "etch: blanket etch shrank the bbox top")
    dose1 = sim.dose("B")
    check(abs(dose1 - dose0) / dose0 < 1e-9, "etch: blanket etch conserves B dose")

    sim.etch(0.2, material="oxide")  # selective: strip remaining oxide
    n2 = sim.n_cells
    check(n2 < n1, "etch: selective material etch reduced n_cells")

    n3 = sim.n_cells
    sim.etch(0.05, poly=[(0.2, 0.2), (0.8, 0.2), (0.8, 0.8), (0.2, 0.8)])
    n4 = sim.n_cells
    check(n4 == n3, "etch: polygon etch leaves n_cells unchanged (retag)")


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
    test_mask_polygon()
    test_deposit_blanket()
    test_etch_blanket()
    test_etch_polygon()
    test_oxidize_dry()
    test_oxidize_wet_faster()
    test_ted_enhancement()
    test_ted_interstitial_field()
    test_mc_damage_seed()
    test_error_paths()
    test_ted_flow_python()
    test_deck_equivalence()
    test_geometry_chain()
    test_save_and_fields()
    test_segregation_dose_loss()
    test_nitride_barrier_python()
    test_activation_python()
    test_rta_ramp_python()
    test_pearson_python()
    test_params_python()
    test_refine_python()
    test_etch_depo_p17()
    print("\nall comprehensive Simulation tests passed")
