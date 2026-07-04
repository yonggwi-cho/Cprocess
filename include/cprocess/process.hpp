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
                          bool seed_damage = false, std::ostream* log = nullptr);

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
double oxidize(SimState& st, double time_s, double temp_k, bool wet = false,
               std::ostream* log = nullptr);

// Adaptively splits mesh edges where `species` has steep concentration
// gradients (M-2). rel_grad_thresh: relative concentration difference across
// a face that triggers refinement (dimensionless, default 0.5 at the
// pybind/Python layer). max_passes: split-pass cap (default 2). Every field
// in st.fields (species included) is carried through the splits by exact
// parent-cell copy; the stack mesh (during photo) is never touched.
void refine(SimState& st, const std::string& species, double rel_grad_thresh,
           int max_passes, std::ostream* log = nullptr);

void save(SimState& st, const std::string& path, std::ostream* log = nullptr);

// Binary CPRC1 state save/load (P1-11). save_state throws if a photoresist
// stack is present (strip it first) -- the stack is never serialized.
// load_state rebuilds mesh topology via finalize() and resets stack members.
void save_state(SimState& st, const std::string& path, std::ostream* log = nullptr);
void load_state(SimState& st, const std::string& path, std::ostream* log = nullptr);

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
