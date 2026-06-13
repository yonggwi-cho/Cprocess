"""
cprocess — 3-D semiconductor process simulator
===============================================

High-level interface::

    import cprocess as cp
    import numpy as np

    # Run a full process deck (SUPREM-style)
    st = cp.SimState()
    log = cp.run_deck('''
        mesh box xmax=0.4um ymax=0.4um zmax=0.8um nx=8 ny=8 nz=40
        init species=P conc=1e15
        implant species=B energy=50keV dose=1e13 method=mc ions=50000
        diffuse time=10min temp=1000C
    ''', st)
    print(log)

    B   = st.get_field('B')    # np.ndarray, shape (n_cells,), [cm^-3]
    xyz = st.mesh.cell_cent    # shape (n_cells, 3), [cm]
    vol = st.mesh.cell_vol     # shape (n_cells,), [cm^3]

Low-level interface::

    mesh = cp.make_box_mesh(0, 0.4*cp.um, 0, 0.4*cp.um, 0, 0.8*cp.um, 8, 8, 40)
    mask = np.ones(mesh.n_cells, dtype=np.uint8)

    p = cp.McImplantParams()
    p.set_dopant('B')
    p.dose = 1e13
    p.energy_kev = 50
    p.ions = 50000
    stats, conc, damage = cp.apply_mc_implant(mesh, mask, p)
    print(stats)            # McImplantStats(deposited=..., Rp=...nm, dRp=...nm)

    opts = cp.DiffuseOpts()
    opts.temp = 1273.15     # 1000 °C in K
    opts.time = 600         # 10 min

    st = cp.SimState()
    # (attach mesh and fields manually or via run_deck)
"""

from ._cprocess import (  # noqa: F401
    # Core objects
    Mesh,
    SimState,

    # Mesh construction
    make_box_mesh,

    # Deck runner
    run_deck,
    run_deck_file,

    # Analytic implant
    ImplantParams,
    apply_implant,

    # MC BCA implant
    McImplantParams,
    McImplantStats,
    apply_mc_implant,

    # Diffusion
    DiffuseOpts,
    DirichletBC,
    diffuse,

    # Output
    save_vtu,
    save_vtu_fields,

    # Unit helpers [in cm or SI-compatible]
    um,
    nm,
    min,
)

# Celsius to Kelvin helper
def celsius(t): return t + 273.15
