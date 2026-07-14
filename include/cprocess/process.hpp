#pragma once
#include <iosfwd>
#include <map>
#include <string>

#include "deck.hpp"
#include "diffusion.hpp"
#include "implant.hpp"
#include "mc_implant.hpp"

namespace cp {

// Native process operations on a SimState. These contain the canonical logic
// for every process step; both the text deck (run_deck) and the Python
// bindings drive them, so the two front-ends always agree. All lengths are in
// cm, energies in keV, times in seconds, temperatures in K — the same
// conventions as the rest of the C++ core. Pass log=nullptr to silence output.
namespace proc {

// Silicon cells (region material "silicon"/"si", or untagged) get mask=1.
std::vector<char> silicon_mask(const SimState& st);

// Resolve a region by name or numeric tag; returns the tag (throws on error).
int resolve_region(const SimState& st, const std::string& name_or_tag);

void mesh_box(SimState& st, double x0, double x1, double y0, double y1,
              double z0, double z1, int nx, int ny, int nz,
              std::ostream* log = nullptr);
void mesh_gmsh(SimState& st, const std::string& file, double scale,
               std::ostream* log = nullptr);

// Set material on a region tag, or on all regions when tag < 0.
void set_region(SimState& st, const std::string& material, int tag = -1,
                std::ostream* log = nullptr);

// Initialize a uniform background concentration (silicon cells only).
void init(SimState& st, const std::string& species, double conc,
          int region = -1, std::ostream* log = nullptr);

// Analytic Gaussian implant. If energy_kev > 0, rp/drp are looked up from the
// dopant range table; otherwise rp/drp must be given. Adds to the field.
double implant_gauss(SimState& st, const std::string& species, double dose,
                     double energy_kev, double rp, double drp, double drl,
                     bool has_window, double x1, double x2, double y1, double y2,
                     bool seed_damage = false,
                     const std::string& profile = "gauss",
                     std::ostream* log = nullptr);

// Monte Carlo (BCA) implant. When a photoresist stack is present (after
// photo()), transport runs through the full physical stack and the profile is
// transferred back onto the working mesh.
// seed_damage=true && channeling=true seeds the "I" field from the MC's own
// Kinchin-Pease damage array x kFrenkelSurvival (capped at
// kAmorphizationDensity); channeling=false falls back to the "+1" model
// (no MC damage is produced without channeling).
McImplantStats implant_mc(SimState& st, const std::string& species, double dose,
                          double energy_kev, long long ions, double tilt_deg,
                          double rotation_deg, unsigned long long seed,
                          int threads, bool channeling, bool has_window,
                          double x1, double x2, double y1, double y2,
                          bool lateral_wrap = false, bool seed_damage = false,
                          std::ostream* log = nullptr);

// Photoresist lithography.
void photo(SimState& st, double thickness, int nz_add = 4,
           std::ostream* log = nullptr);
// Rectangular mask opening (convenience wrapper around mask_polygon).
void mask(SimState& st, double x1, double x2, double y1, double y2,
          std::ostream* log = nullptr);
// Arbitrary 2-D polygon mask opening. `poly` is a list of (x,y) vertices in
// cm defining a closed polygon; winding-number point-in-polygon test is used.
// Resist cells whose centroid projects inside the polygon are opened (vacuum).
void mask_polygon(SimState& st,
                  const std::vector<std::pair<double,double>>& poly,
                  std::ostream* log = nullptr);
void strip(SimState& st, std::ostream* log = nullptr);

// Deposit a blanket or polygon-shaped film of `material` (oxide/nitride/poly)
// on the current top surface. `poly` restricts deposition to cells whose XY
// centroid is inside the polygon; empty vector means blanket (full surface).
// `thickness` in cm, returned mesh is the new top after deposition.
void deposit(SimState& st, const std::string& material,
             double thickness, int nz_add = 2,
             const std::vector<std::pair<double,double>>& poly = {},
             std::ostream* log = nullptr);

// Etch the top surface down by `depth` (cm).
// poly empty (or < 3 vertices): BLANKET etch — cells with centroid above
// (bbox.hi.z - depth) and matching `material` are physically REMOVED from
// the mesh (node compaction + identity field transfer), shrinking the mesh
// and bbox. poly given: POLYGON etch — matching cells inside the polygon
// keep their old behavior of being retagged to "gas" with zeroed fields (the
// mesh stays box-shaped; true removal of an interior column would produce a
// non-manifold surface that break extend_mesh_exact()/infer_box_dims() on
// the next deposit — that generalization is P2-5's scope).
// material: empty means all non-gas materials; otherwise only cells whose
// region_material matches (case-insensitive, "si"/"silicon" aliased) are
// affected — other cells in the depth band are left untouched (selectivity).
// Throws if a photoresist stack is present (strip first).
void etch(SimState& st, double depth,
          const std::vector<std::pair<double,double>>& poly = {},
          const std::string& material = "",
          std::ostream* log = nullptr);

// P2-5: rate/time level-set etch. Uses the cp::levelset toolkit (levelset.hpp)
// to advect a signed-distance field from the exposed surface, so (unlike the
// geometric etch() above) an isotropic etch produces a genuine lateral
// undercut under a mask overhang, and a vertical/anisotropic etch removes
// material only where the local surface normal faces "up" (away from
// remaining solid, into gas). `rates`: material name (case-insensitive,
// "si"/"silicon" aliased) -> etch rate in cm/s; a material absent from the
// map has rate 0 (untouched). `poly` (cm, closed, empty = blanket) further
// restricts which (x,y) columns are exposed to the given rates -- material
// outside the polygon is not attacked even if its rate is nonzero. This
// coexists with etch() (geometric depth/poly) -- neither replaces the other.
void etch_rate(SimState& st, const std::map<std::string, double>& rates,
               double time_s, bool isotropic = true,
               const std::vector<std::pair<double,double>>& poly = {},
               std::ostream* log = nullptr);

// P2-5: conformal (isotropic level-set) deposit of `material`, `thickness_cm`
// thick, measured normal-to-surface everywhere -- including down sidewalls of
// any existing step/trench, unlike deposit()'s purely vertical film growth.
void deposit_conformal(SimState& st, const std::string& material,
                       double thickness_cm, std::ostream* log = nullptr);

// P3-b: epitaxial growth of `thickness_cm` of silicon on the exposed silicon
// top surface, with in-situ uniform doping `doping` (species symbol -> cm^-3)
// in the newly grown cells only. Geometry reuses deposit()'s extend_mesh_exact
// + retag + layer_stack mechanism with material "silicon". When anneal ==
// true (default) the growth thermal budget is applied by one automatic
// proc::diffuse_ted(temp_k, time_s) call after growth, so substrate dopants
// back-diffuse into the epi layer. No growth-rate model: thickness and time
// are both caller-given (P3_overview: confirmed policy).
void epitaxy(SimState& st, double thickness_cm, double temp_k,
             double time_s,
             const std::map<std::string, double>& doping = {},
             bool anneal = true, std::ostream* log = nullptr);

// Dirichlet boundary conditions for diffusion.
void add_bc(SimState& st, const std::string& species, int patch, double conc,
            std::ostream* log = nullptr);
void clear_bc(SimState& st, std::ostream* log = nullptr);

void diffuse(SimState& st, const DiffuseOpts& opts, std::ostream* log = nullptr);

// Transient enhanced diffusion. Uses the interstitial excess accumulated in the
// "I" field (seeded by implants when damage=true, per the "+1" model) to
// enhance dopant diffusivity, capturing the initial fast-diffusion transient.
// The interstitial field decays via bulk recombination and the surface sink.
void diffuse_ted(SimState& st, const DiffuseOpts& opts,
                 std::ostream* log = nullptr);

// Blanket 1D vertical thermal oxidation of the exposed top surface.
// time_s in seconds, temp_k in K. Grows/extends the SiO2 layer above the
// silicon per Deal-Grove; consumes 0.44*dx_ox of Si and raises the outer
// surface by 0.56*dx_ox. Returns the new total oxide thickness in cm.
// P2-3 (OED): when ParamDB "oed.theta" > 0 (default 0.01), the call
// internally sub-steps into N=10 (grow -> inject excess "I" at the growing
// interface -> relax via diffuse_ted) increments, so a single oxidize()
// call already exhibits oxidation-enhanced diffusion for any dopant field
// present. "oed.theta" == 0 disables this and reproduces the original
// (P1-6) geometry-only behavior exactly (no diffusion at all inside
// oxidize()).
double oxidize(SimState& st, double time_s, double temp_k, bool wet = false,
               std::ostream* log = nullptr);

// 2D/3D LOCOS-style oxidation (P2-4): lateral bird's-beak encroachment under
// a nitride mask, via a genuine steady-state oxidant-diffusion solve
// (cp::solve_oxidant) each sub-step instead of a single blanket Deal-Grove
// number. Requires a nitride-masked region to exist in the mesh (region
// material "nitride" somewhere on the current top surface) -- for a fully
// blanket (no mask) oxidation, use oxidize() instead, which is unchanged and
// remains the fast 1D path.
//
// API design note: this is a separate entry point rather than a branch
// inside oxidize() itself. oxidize() is a heavily-tuned, sub-stepped
// incremental-retag implementation (P1-6/P2-3) whose numerics (deferred
// geometry realization, OED "I" injection caps, diffuse_ted solver
// tolerances) are calibrated against its own test suite; adding a second,
// structurally different geometry model (column-wise, solver-driven
// interface heights) as an in-place branch would risk destabilizing that
// tuning for zero benefit to the blanket case. Keeping the 2D path as its
// own function guarantees the blanket 1D flow can never regress.
//
// Implementation model: the mesh must be a make_box_mesh-style regular grid
// (as built/extended by mesh_box/deposit/etch/oxidize elsewhere in this
// API). Growth is resolved per (x,y) grid column: each sub-step, solve_oxidant
// is run over the oxide+nitride sub-mesh (nitride cells get a small but
// nonzero oxidant diffusivity, `ox2d.nitride_leak` (default 1e-3) x the
// oxide value -- this is what lets a little oxidant bleed laterally under
// the mask edge, combined with lateral diffusion from the open field, to
// produce the bird's-beak taper), giving each interface face a local growth
// velocity; these are averaged per column to get a column-local dx_ox, then
// realized as a per-column interface height (0.44*dx_ox consumed into Si,
// 0.56*dx_ox raising the outer surface) via the same extend+retag+
// repair_quality machinery oxidize() uses, generalized from a single
// interface height to a per-column height field -- the box mesh's regular
// structure makes this equivalent to genuine per-node normal motion (each
// node belongs to exactly one column) without needing the general
// unstructured ALE traversal. If no oxide exists anywhere yet, a thin
// uniform "native" seed layer (`ox2d.seed_ox_um` default 0.001 um) is grown
// over the open (non-nitride) columns first so solve_oxidant has cells to
// operate on.
double oxidize_2d(SimState& st, double time_s, double temp_k, bool wet = false,
                  std::ostream* log = nullptr);

// Blanket silicidation (P3-d): converts a previously deposited blanket
// metal film ("nickel" or "titanium") on the exposed Si top surface into
// its silicide ("nisi" / "tisi2") by diffusion-limited growth
//   x^2 = x0^2 + B(T)*time_s,  B = b0*exp(-eb/kT)   (ParamDB, per phase).
// Volume bookkeeping (all retag-only, no mesh rebuild -- the surface
// *recedes* by (rsi+rmet-1)*dx): growing dx of silicide consumes
// rsi*dx of Si (interface moves down) and rmet*dx of metal; the excess
// band at the top is retagged "gas". Growth stops when the metal is
// exhausted (x capped at x0 + t_metal/rmet). Returns the new total
// silicide thickness in cm. metal: "nickel"|"ni"|"titanium"|"ti".
double silicide(SimState& st, const std::string& metal, double temp_k,
                double time_s, std::ostream* log = nullptr);

// Linear-elastic FEM mechanics solve (P2-6): builds a per-cell eigenstrain
// load from thermal mismatch (alpha_m * dT, dT = temp_k - 300K) and intrinsic
// film stress (ParamDB "mech.sigma0.<material>", isotropic), applies
// Dirichlet BCs (roller only at the three "min" faces -- zmin/xmin/ymin,
// each pinning just its own normal displacement component -- which removes
// all 6 rigid-body modes while staying compatible with pure thermal
// expansion; pinning both sides of an axis, or all 3 components at zmin,
// would suppress thermal expansion outright), solves K u = f via cg_ilu0,
// computes per-cell Voigt stress
// sigma = D*(B*u - eps0), then applies one dt_s-worth of Maxwell relaxation
// (cp::maxwell_update, using ParamDB "mech.tau.<material>" in seconds; 0 or
// absent means purely elastic/no relaxation). Results are written to
// st.fields["sxx"], "syy", "szz", "sxy", "syz", "sxz" (dyn/cm^2, the CGS
// stress unit consistent with the rest of the C++ core).
// P3-f: if st.fields["Ge"] is present, its composition (x_Ge = C_Ge/kNSi)
// contributes an additional Vegard's-law eigenstrain (ParamDB
// "sige.eps0_coef", default 0.042) added in parallel with the thermal/
// intrinsic terms above; this is the strain *source* only. To have that
// strain actually perturb diffusion, call set_param("stress.couple", 1)
// before diffuse()/diffuse_ted() so P3-e's DiffuseOpts::pressure hydrostatic
// coupling picks up the resulting stress field from this call.
void mechanics(SimState& st, double temp_k, double dt_s,
              std::ostream* log = nullptr);

// Adaptively splits mesh edges where `species` has steep concentration
// gradients (M-2). rel_grad_thresh: relative concentration difference across
// a face that triggers refinement (dimensionless, default 0.5 at the
// pybind/Python layer). max_passes: split-pass cap (default 2). Every field
// in st.fields (species included) is carried through the splits by exact
// parent-cell copy; the stack mesh (during photo) is never touched.
// axis (M-6): "" (default) keeps the M-2 gradient-driven behavior exactly.
// "x"|"y"|"z" switches to direction-aligned refinement (refine_anisotropic):
// the gradient indicator is not used, `direction` is the corresponding unit
// axis, and `rel_grad_thresh` is reused as align_thresh (0.5 is a sane
// default for both). `species` is only used to confirm the field exists (for
// field-transfer bookkeeping), not to drive the indicator. Any other axis
// value throws std::runtime_error.
void refine(SimState& st, const std::string& species, double rel_grad_thresh,
           int max_passes, const std::string& axis = "",
           std::ostream* log = nullptr);

void save(SimState& st, const std::string& path, std::ostream* log = nullptr);

// Binary CPRC1 state save/load (P1-11). save_state throws if a photoresist
// stack is present (strip it first) -- the stack is never serialized.
// load_state rebuilds mesh topology via finalize() and resets stack members.
void save_state(SimState& st, const std::string& path, std::ostream* log = nullptr);
void load_state(SimState& st, const std::string& path, std::ostream* log = nullptr);

// P3-h: device-simulator export. Writes <prefix>.vtu (cell data as
// proc::save + node-averaged dopant/active/NetDoping point data) and
// <prefix>.meta.json (region/material table, boundary patch names, unit
// system, species list, ND-NA sign convention). Throws if a photoresist
// stack is present (strip first) or on I/O error.
void export_device(SimState& st, const std::string& path_prefix,
                   std::ostream* log = nullptr);

// 1-D depth profile of `species` at column (x_cm, y_cm): every non-gas cell
// whose xy-bbox contains the point, sorted by cell_cent.z ascending. Throws
// if the species field does not exist; returns an empty vector if no cell's
// xy-bbox contains the point.
std::vector<std::pair<double, double>> profile1d(const SimState& st,
                                                  const std::string& species,
                                                  double x_cm, double y_cm);

// Per-cell electrically active concentration of `species` at temp_k [K]
// (solid-solubility clamp, P1-3). temp_k <= 0 selects st.last_temp. Throws if
// the field does not exist.
std::vector<double> active_field(const SimState& st, const std::string& species,
                                 double temp_k = -1.0, std::ostream* log = nullptr);

// GDSII layout input (P2-7). Reads BOUNDARY polygons on `layer` (layer < 0
// means all layers) from a minimal dependency-free GDSII stream reader
// (cp::read_gds); returned polygons are in cm, ready to pass to
// mask_polygon()/deposit(poly=)/etch(poly=). Throws on missing/malformed
// files (see cp::read_gds).
std::vector<std::vector<std::pair<double, double>>> load_gds(
    const std::string& path, int layer, std::ostream* log = nullptr);

// ── Runtime parameter overrides (P1-10, cp::ParamDB) ──
// Override any physical parameter consumed by cp::materials.cpp / the TED
// loop / seed_interstitials at call time (raw core units). No key validation
// -- an unknown key is silently harmless. `st` is unused (kept for API
// consistency with the other proc:: functions).
void set_param(SimState& st, const std::string& key, double value,
               std::ostream* log = nullptr);
// Returns the override for `key` if set, else `fallback`.
double get_param(const std::string& key, double fallback = 0.0);
// Currently-set overrides (not the full supported key set).
std::map<std::string, double> list_params();

}  // namespace proc
}  // namespace cp
