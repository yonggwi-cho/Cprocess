"""
cprocess — 3-D semiconductor process simulator
===============================================

Native Python interface (recommended)::

    import cprocess as cp

    sim = cp.Simulation()
    sim.mesh(x=1.0, y=1.0, z=0.5, nx=8, ny=8, nz=4)   # micrometres
    sim.region("silicon")
    sim.init("B", 1e15)

    sim.implant("B", dose=5e12, energy=40, mc=True)    # channel-stop

    sim.photo(resist=0.4)                              # physical resist
    sim.mask(x1=0.0,  x2=0.35)                          # open source
    sim.mask(x1=0.65, x2=1.0)                           # open drain
    sim.implant("P", dose=5e15, energy=30, mc=True)
    sim.strip()

    sim.diffuse(time=30, temp=1000)                    # min, Celsius
    sim.save("nmos.vtu")

    P   = sim.field("P")          # np.ndarray [cm^-3]
    xyz = sim.cell_centroids      # np.ndarray [n_cells, 3], micrometres

Units: micrometres, keV, minutes, Celsius. Methods chain (each returns the
Simulation), except implant() which returns the implant result.

Text deck (legacy, still supported)::

    st = cp.SimState()
    cp.run_deck('''
        mesh box xmax=0.4um ymax=0.4um zmax=0.8um nx=8 ny=8 nz=40
        init species=P conc=1e15
        implant species=B energy=50keV dose=1e13 method=mc ions=50000
        diffuse time=10min temp=1000C
    ''', st)

Low-level bindings (make_box_mesh, apply_mc_implant, diffuse, ...) remain
available for custom workflows.
"""

from .simulation import Simulation, ImplantResult, UM, NM, MIN  # noqa: F401

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
