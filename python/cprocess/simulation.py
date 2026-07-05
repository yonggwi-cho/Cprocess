"""
Native Python process-simulation interface.

`Simulation` is the recommended front-end: it drives the C++ process engine
directly, with friendly engineering units (micrometres, keV, minutes, Celsius)
instead of the text deck. Every method maps one-to-one onto a process step.

Example — masked NMOS source/drain through a physical photoresist gate::

    import cprocess as cp

    sim = cp.Simulation()
    sim.mesh(x=1.0, y=1.0, z=0.5, nx=8, ny=8, nz=4)   # micrometres
    sim.region("silicon")
    sim.init("B", 1e15)                                # p-type substrate

    sim.implant("B", dose=5e12, energy=40, mc=True)    # channel-stop

    sim.photo(resist=0.4)                              # 0.4 µm resist
    sim.mask(x1=0.0,  x2=0.35)                          # open source
    sim.mask(x1=0.65, x2=1.0)                           # open drain
    sim.implant("P", dose=5e15, energy=30, mc=True, ions=80000)
    sim.strip()

    sim.diffuse(time=30, temp=1000)                    # min, Celsius
    sim.save("nmos.vtu")

    P = sim.field("P")          # numpy array [cm^-3]
    xyz = sim.cell_centroids    # numpy [n_cells, 3], micrometres
"""

from __future__ import annotations

import numpy as np

from . import _cprocess as _c

UM = 1e-4   # micrometre in cm
NM = 1e-7   # nanometre in cm
MIN = 60.0  # minute in seconds


def _celsius_to_k(t_c: float) -> float:
    return t_c + 273.15


class ImplantResult:
    """Outcome of a Monte-Carlo implant step."""

    def __init__(self, stats, log: str):
        self.stats = stats
        self.log = log

    @property
    def deposited(self) -> int:
        return self.stats.deposited

    @property
    def rp_nm(self) -> float:
        return self.stats.rp * 1e7

    @property
    def drp_nm(self) -> float:
        return self.stats.drp * 1e7

    def __repr__(self) -> str:
        return (f"ImplantResult(deposited={self.deposited}, "
                f"Rp={self.rp_nm:.1f} nm, dRp={self.drp_nm:.1f} nm)")


class Simulation:
    """A 3-D process simulation driven natively from Python.

    Lengths are micrometres, energies keV, times minutes, temperatures Celsius
    unless noted. The underlying C++ engine works in cm/keV/s/K; conversions are
    handled here.
    """

    def __init__(self, verbose: bool = False):
        self._st = _c.SimState()
        self.verbose = verbose
        self.log = ""

    # -- internal --------------------------------------------------------------
    def _emit(self, text: str) -> None:
        if text:
            self.log += text
            if self.verbose:
                print(text, end="")

    # -- geometry --------------------------------------------------------------
    def mesh(self, x: float, y: float, z: float,
             nx: int, ny: int, nz: int,
             x0: float = 0.0, y0: float = 0.0, z0: float = 0.0) -> "Simulation":
        """Create a box mesh. Extents and origin in micrometres."""
        self._emit(_c.proc_mesh_box(self._st,
                   x0 * UM, x * UM, y0 * UM, y * UM, z0 * UM, z * UM,
                   int(nx), int(ny), int(nz)))
        return self

    def mesh_gmsh(self, path: str, scale_um: float = 1.0) -> "Simulation":
        """Load a Gmsh mesh; `scale_um` converts file units to micrometres."""
        self._emit(_c.proc_mesh_gmsh(self._st, path, scale_um * UM))
        return self

    def region(self, material: str, tag: int = -1) -> "Simulation":
        """Set a region's material (tag<0 = all regions).

        material ∈ {silicon, oxide, nitride, poly, gas}. Only silicon cells are
        implanted and diffused.
        """
        self._emit(_c.proc_set_region(self._st, material, int(tag)))
        return self

    # -- doping ----------------------------------------------------------------
    def init(self, species: str, conc: float, region=None) -> "Simulation":
        """Uniform background doping [cm^-3] (silicon cells)."""
        reg = -1 if region is None else (
            region if isinstance(region, int)
            else _c.proc_resolve_region(self._st, str(region)))
        self._emit(_c.proc_init(self._st, species, float(conc), int(reg)))
        return self

    def implant(self, species: str, dose: float,
                energy: float = 0.0, *,
                mc: bool = False,
                rp: float = 0.0, drp: float = 0.0, drl: float = 0.0,
                ions: int = 100000, tilt: float = 0.0, rotation: float = 0.0,
                seed: int = 1, threads: int = 0, channeling: bool = False,
                window=None, damage: bool = False, profile: str = "gauss"):
        """Ion implant.

        dose [cm^-2], energy [keV]. `mc=True` selects Monte-Carlo BCA, otherwise
        an analytic Gaussian (give energy, or rp/drp/drl in micrometres).
        `window=(x1,x2,y1,y2)` in micrometres restricts a geometric mask; for a
        *physical* resist mask use photo()/mask() instead. If a resist stack is
        present, an MC implant automatically transports through it.
        `damage=True` seeds excess self-interstitials into the "I" field, which
        diffuse(ted=True) then uses for transient enhanced diffusion. For
        `mc=True` with `channeling=True`, the seed comes from the MC's own
        Kinchin-Pease damage field (x kFrenkelSurvival, capped at
        kAmorphizationDensity); otherwise (Gaussian, or MC without channeling)
        the "+1" model is used (dopant profile copied into "I").
        `profile="pearson"` selects a Pearson-IV depth profile (analytic
        implant only; requires `energy=`, not rp/drp) built from the 4-moment
        table (Rp, dRp, gamma, beta); it falls back to Gaussian when the
        moments don't satisfy the Type-IV validity condition.
        """
        if mc and profile != "gauss":
            raise ValueError("profile applies to the analytic implant only")
        has_window = window is not None
        x1, x2, y1, y2 = (window if has_window else (0, 0, 0, 0))
        if mc:
            stats, log = _c.proc_implant_mc(self._st, species, float(dose),
                float(energy), int(ions), float(tilt), float(rotation),
                int(seed), int(threads), bool(channeling), has_window,
                x1 * UM, x2 * UM, y1 * UM, y2 * UM, bool(damage))
            self._emit(log)
            return ImplantResult(stats, log)
        atoms, log = _c.proc_implant_gauss(self._st, species, float(dose),
            float(energy), rp * UM, drp * UM, drl * UM, has_window,
            x1 * UM, x2 * UM, y1 * UM, y2 * UM, bool(damage), profile)
        self._emit(log)
        return atoms

    # -- lithography -----------------------------------------------------------
    def photo(self, resist: float, nz: int = 4) -> "Simulation":
        """Deposit a blanket photoresist layer `resist` micrometres thick."""
        self._emit(_c.proc_photo(self._st, resist * UM, int(nz)))
        return self

    def mask(self, x1: float, x2: float,
             y1: float = None, y2: float = None) -> "Simulation":
        """Expose and develop a rectangular window in the resist (micrometres).

        y-range defaults to the full device width.
        """
        bb = self.bbox()
        yy1 = bb[0][1] if y1 is None else y1 * UM
        yy2 = bb[1][1] if y2 is None else y2 * UM
        self._emit(_c.proc_mask(self._st, x1 * UM, x2 * UM, yy1, yy2))
        return self

    def mask_polygon(self, poly) -> "Simulation":
        """Expose a polygon-shaped window in the resist.

        poly: sequence of (x, y) tuples in micrometres, e.g.::

            sim.mask_polygon([(0.1, 0.1), (0.3, 0.0), (0.5, 0.2), (0.2, 0.4)])

        The polygon is automatically closed (last vertex connects to first).
        """
        poly_cm = [(x * UM, y * UM) for x, y in poly]
        self._emit(_c.proc_mask_polygon(self._st, poly_cm))
        return self

    def load_gds(self, path: str, layer: int = -1, scale: float = 1.0):
        """Read BOUNDARY polygons from a GDSII stream file.

        path: path to a .gds file.
        layer: GDS layer number to keep (-1: all layers, default).
        scale: extra multiplier applied on top of the GDS UNITS-derived
            scale (default 1.0 -- normally leave this alone).

        Returns a list of polygons, each a list of (x, y) tuples in
        micrometres, suitable to pass directly to mask_polygon()/
        deposit(poly=...)/etch(poly=...). This is a query method: it
        returns data, not self.
        """
        polys_cm = _c.load_gds(path, int(layer))
        return [[(x / UM * scale, y / UM * scale) for x, y in poly]
                for poly in polys_cm]

    def deposit(self, material: str, thickness: float,
                nz: int = 2, poly=None) -> "Simulation":
        """Deposit a film on the top surface. Multi-layer deposits are

        supported (e.g. deposit('oxide', ...).deposit('nitride', ...)): each
        call pushes a new entry onto the internal layer stack, most recent
        first, so earlier films keep their own material identity.

        material: 'oxide', 'nitride', 'poly', 'silicon', ...
        thickness: film thickness in micrometres.
        poly: optional sequence of (x, y) tuples in micrometres that restricts
              the deposit to a polygon footprint; omit for blanket deposition.
        """
        poly_cm = [(x * UM, y * UM) for x, y in poly] if poly else []
        self._emit(_c.proc_deposit(self._st, material, thickness * UM,
                                   int(nz), poly_cm))
        return self

    def etch(self, depth: float, poly=None, material: str = "") -> "Simulation":
        """Etch the top surface down by `depth` micrometres.

        poly: optional sequence of (x, y) tuples in micrometres that restricts
              the etch to a polygon footprint; omit for blanket etch.
        material: optional material name (e.g. 'oxide') to etch selectively;
              empty (default) etches all non-gas materials.

        A blanket etch (poly omitted) physically REMOVES the etched cells
        from the mesh (node compaction + identity field transfer for the
        surviving cells) -- `n_cells` and the bbox shrink accordingly. A
        polygon etch keeps the legacy behavior of retagging matching cells
        as 'gas' with zeroed concentrations, leaving the mesh topology (and
        n_cells) unchanged -- true removal of an interior column would
        produce a non-manifold top surface that breaks the box-mesh
        assumptions of the next deposit().
        """
        poly_cm = [(x * UM, y * UM) for x, y in poly] if poly else []
        self._emit(_c.proc_etch(self._st, depth * UM, poly_cm, material))
        return self

    def etch_rate(self, rates: dict, time: float, isotropic: bool = True,
                 poly=None) -> "Simulation":
        """Level-set rate/time etch (P2-5).

        rates: {material_name: rate_um_per_min}, e.g. {"silicon": 0.1}.
              Materials absent from the dict are not attacked (rate 0).
        time: etch duration in minutes.
        isotropic: True undercuts laterally under mask overhangs (a
              standalone signed-distance-field advection, unlike the
              vertical-only geometric etch()); False removes material only
              where the local surface normal faces away from remaining
              solid (a vertical/RIE approximation).
        poly: optional sequence of (x, y) tuples in micrometres restricting
              which (x, y) columns are exposed to the given rates (empty =
              blanket).

        Coexists with etch() (geometric depth/poly); neither replaces the
        other.
        """
        rates_cm_s = {mat: r * UM / MIN for mat, r in rates.items()}
        poly_cm = [(x * UM, y * UM) for x, y in poly] if poly else []
        self._emit(_c.proc_etch_rate(self._st, rates_cm_s, time * MIN,
                                      bool(isotropic), poly_cm))
        return self

    def deposit_conformal(self, material: str, thickness: float) -> "Simulation":
        """Conformal (isotropic level-set) deposit of `material`.

        thickness: film thickness in micrometres, measured normal-to-
              surface everywhere -- including down the sidewalls of an
              existing step/trench, unlike deposit()'s purely vertical
              film growth.

        Coexists with deposit(); neither replaces the other.
        """
        self._emit(_c.proc_deposit_conformal(self._st, material,
                                              thickness * UM))
        return self

    def oxidize(self, time: float, temp: float, *, wet: bool = False) -> "Simulation":
        """Blanket thermal oxidation of the exposed silicon top surface.

        time in minutes, temp in Celsius. Grows SiO2 per Deal-Grove
        (<100> Si); the surface rises by 0.56x and silicon is consumed
        by 0.44x of the grown oxide thickness.
        """
        self._emit(_c.proc_oxidize(self._st, time * MIN,
                                   _celsius_to_k(temp), bool(wet)))
        return self

    def oxidize_2d(self, time: float, temp: float, *, wet: bool = False) -> "Simulation":
        """2D/3D LOCOS oxidation (P2-4): bird's-beak lateral encroachment.

        time in minutes, temp in Celsius. Requires a nitride-masked region
        to already be present in the mesh (deposit("nitride", ...) then
        etch an opening); raises RuntimeError otherwise -- use oxidize()
        for blanket (1D) oxidation. Solves a steady-state oxidant-diffusion
        model each sub-step instead of a single blanket Deal-Grove number,
        so the oxide tapers laterally under the mask edge.
        """
        self._emit(_c.proc_oxidize_2d(self._st, time * MIN,
                                      _celsius_to_k(temp), bool(wet)))
        return self

    def mechanics(self, temp: float, time: float) -> "Simulation":
        """Linear-elastic FEM mechanics solve (P2-6).

        temp in Celsius, time in minutes. Builds a per-cell eigenstrain load
        from thermal mismatch (dT = temp - 300K reference) and any intrinsic
        film stress (ParamDB "mech.sigma0.<material>"), applies Dirichlet BCs
        (zmin fixed, lateral faces roller), solves for the displacement field
        via cg_ilu0, and writes the per-cell Voigt stress (dyn/cm^2) to fields
        "sxx", "syy", "szz", "sxy", "syz", "sxz". Then applies `time` worth of
        Maxwell relaxation (ParamDB "mech.tau.<material>" in seconds; 0/unset
        means purely elastic).
        """
        self._emit(_c.proc_mechanics(self._st, _celsius_to_k(temp), time * MIN))
        return self

    def strip(self) -> "Simulation":
        """Remove all remaining photoresist."""
        self._emit(_c.proc_strip(self._st))
        return self

    def refine(self, species: str, threshold: float = 0.5,
              passes: int = 2) -> "Simulation":
        """Adaptively split mesh edges where `species` has steep gradients.

        `threshold` is a dimensionless relative concentration-difference
        threshold across a face (no unit conversion). `passes` caps the
        number of split rounds.
        """
        self._emit(_c.proc_refine(self._st, species, float(threshold),
                                  int(passes)))
        return self

    # -- diffusion -------------------------------------------------------------
    def bc(self, species: str, patch: str, conc: float) -> "Simulation":
        """Fixed-concentration Dirichlet BC on a named boundary patch."""
        p = _c.find_patch(self._st, patch)
        if p < 0:
            raise ValueError(f"unknown patch '{patch}'")
        self._emit(_c.proc_add_bc(self._st, species, p, float(conc)))
        return self

    def clear_bc(self) -> "Simulation":
        _c.proc_clear_bc(self._st)
        return self

    def diffuse(self, time: float, temp: float = None, *,
                dt: float = 0.0, field_enh: bool = True,
                nonortho: bool = True, ted: bool = False,
                activation: bool = True, ramp=None) -> "Simulation":
        """Anneal: `time` in minutes, `temp` in Celsius, `dt` in minutes.

        `ted=True` enables transient enhanced diffusion, coupling the excess
        self-interstitials seeded by implant(damage=True) into the dopant
        diffusivity. The enhancement decays as interstitials reach the surface
        sink and recombine, reproducing the initial fast-diffusion transient.

        `activation=True` (default) clamps the concentration entering charge
        neutrality (and hence field enhancement) at the solid solubility
        C_ss(T); set False to treat the full concentration as active.

        `ramp`, if given, is a piecewise-linear RTA temperature profile
        `[(t_min, T_celsius), ...]` with `ramp[0][0] == 0`. When `ramp` is
        given, `temp` is ignored (the profile's first point sets the initial
        temperature); one of `temp` or `ramp` must be given.
        """
        opts = _c.DiffuseOpts()
        opts.time = time * MIN
        if ramp is not None:
            if ramp[0][0] != 0:
                raise ValueError("ramp must start at t=0")
            opts.temp_profile = [(t * MIN, _celsius_to_k(T)) for t, T in ramp]
            opts.temp = _celsius_to_k(ramp[0][1])
        elif temp is not None:
            opts.temp = _celsius_to_k(temp)
        else:
            raise ValueError("give temp or ramp")
        opts.dt = dt * MIN
        opts.field_enh = field_enh
        opts.nonortho = nonortho
        opts.verbosity = 1 if self.verbose else 0
        opts.activation = activation
        if ted:
            self._emit(_c.proc_diffuse_ted(self._st, opts))
        else:
            self._emit(_c.proc_diffuse(self._st, opts))
        return self

    # -- output ----------------------------------------------------------------
    def save(self, path: str) -> "Simulation":
        """Write the mesh and all fields to a ParaView .vtu file."""
        self._emit(_c.proc_save(self._st, path))
        return self

    def save_state(self, path: str) -> "Simulation":
        """Save the full simulation state to a binary CPRC1 file."""
        self._emit(_c.proc_save_state(self._st, path))
        return self

    def load_state(self, path: str) -> "Simulation":
        """Restore a state previously written by save_state()."""
        self._emit(_c.proc_load_state(self._st, path))
        return self

    def profile(self, species: str, x: float, y: float):
        """Depth profile at column (x, y) [um].
        Returns (z_um, conc) numpy arrays sorted by depth."""
        pairs = _c.proc_profile1d(self._st, species, x * UM, y * UM)
        z = np.array([p[0] for p in pairs]) / UM
        c = np.array([p[1] for p in pairs])
        return z, c

    # -- parameter overrides (P1-10) --------------------------------------------
    def set_param(self, key: str, value: float) -> "Simulation":
        """Override a physical parameter (raw core units: cm^2/s, eV, cm^-3,
        1/s, dimensionless -- no unit conversion). e.g.
        set_param("B.d0", 0.074). See materials.hpp for the supported key set;
        unknown keys are silently inert."""
        self._emit(_c.proc_set_param(self._st, key, float(value)))
        return self

    def get_param(self, key: str, fallback: float = 0.0) -> float:
        """Current override for `key`, or `fallback` if unset (raw core units)."""
        return _c.proc_get_param(key, float(fallback))

    # -- accessors -------------------------------------------------------------
    @property
    def state(self):
        return self._st

    @property
    def mesh_obj(self):
        return self._st.mesh

    @property
    def n_cells(self) -> int:
        return self._st.mesh.n_cells

    @property
    def cell_volumes(self) -> np.ndarray:
        """Per-cell volume [cm^3]."""
        return self._st.mesh.cell_vol

    @property
    def cell_centroids(self) -> np.ndarray:
        """Per-cell centroid in micrometres, shape (n_cells, 3)."""
        return self._st.mesh.cell_cent / UM

    def field_names(self):
        return self._st.field_names()

    def field(self, species: str) -> np.ndarray:
        """Concentration field [cm^-3], shape (n_cells,)."""
        return self._st.get_field(species)

    def active(self, species: str, temp: float = None) -> np.ndarray:
        """Electrically active concentration [cm^-3] (solid-solubility clamp).

        temp in Celsius; None uses the temperature of the last diffuse step.
        """
        t_k = -1.0 if temp is None else _celsius_to_k(temp)
        return _c.proc_active_field(self._st, species, t_k)

    def set_field(self, species: str, conc: np.ndarray) -> "Simulation":
        self._st.set_field(species, np.asarray(conc, dtype=float))
        return self

    def dose(self, species: str) -> float:
        """Integrated areal dose in silicon [cm^-2]."""
        c = self.field(species)
        v = self.cell_volumes
        bb = self.bbox()
        area = (bb[1][0] - bb[0][0]) * (bb[1][1] - bb[0][1])
        return float(np.sum(c * v) / area) if area > 0 else 0.0

    def bbox(self):
        """((x0,y0,z0),(x1,y1,z1)) in cm."""
        xyz = self._st.mesh.nodes
        return (xyz.min(axis=0), xyz.max(axis=0))
