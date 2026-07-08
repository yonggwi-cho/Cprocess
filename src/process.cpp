#include "cprocess/process.hpp"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <map>
#include <ostream>
#include <queue>
#include <stdexcept>
#include <utility>

#include "cprocess/ale_mover.hpp"
#include "cprocess/device_export.hpp"
#include "cprocess/fem.hpp"
#include "cprocess/field_transfer.hpp"
#include "cprocess/gds_reader.hpp"
#include "cprocess/gmsh_reader.hpp"
#include "cprocess/levelset.hpp"
#include "cprocess/mechanics.hpp"
#include "cprocess/oxidant_solver.hpp"
#include "cprocess/oxidation.hpp"
#include "cprocess/param_db.hpp"
#include "cprocess/remesh.hpp"
#include "cprocess/state_io.hpp"
#include "cprocess/topology.hpp"
#include "cprocess/vtk_writer.hpp"

namespace cp {
namespace proc {

namespace {

std::string lower(std::string s) {
  std::transform(s.begin(), s.end(), s.begin(),
                 [](unsigned char c) { return std::tolower(c); });
  return s;
}

bool is_silicon(const std::string& mat) {
  const std::string m = lower(mat);
  return m == "silicon" || m == "si";
}

std::string fmt(const char* f, double v) {
  char buf[64];
  std::snprintf(buf, sizeof(buf), f, v);
  return buf;
}

// Build a box mesh covering the same (x,y) footprint as `base` but taller by
// `thickness`, with `nz_add` extra layers stacked above the substrate.
// Infer the (nx, ny, nz) cell counts of a regular box mesh from its nodes.
// For a box mesh built with make_box_mesh(x0,x1,y0,y1,z0,z1,nx,ny,nz) there
// are (nx+1)*(ny+1)*(nz+1) nodes placed on a regular grid; we recover nx/ny/nz
// by counting unique coordinate values on each axis.
static std::tuple<int,int,int> infer_box_dims(const Mesh& m) {
  if (m.nodes.empty()) return {1, 1, 1};
  std::vector<double> xs, ys, zs;
  xs.reserve(m.nodes.size()); ys.reserve(m.nodes.size()); zs.reserve(m.nodes.size());
  for (const Vec3& n : m.nodes) {
    xs.push_back(n.x); ys.push_back(n.y); zs.push_back(n.z);
  }
  // Count unique values with tolerance.
  auto count_unique = [](std::vector<double>& v) {
    std::sort(v.begin(), v.end());
    int cnt = 1;
    for (std::size_t i = 1; i < v.size(); ++i)
      if (v[i] - v[i-1] > 1e-12 * (v.back() - v.front() + 1e-20)) ++cnt;
    return cnt;
  };
  const int nx = std::max(1, count_unique(xs) - 1);
  const int ny = std::max(1, count_unique(ys) - 1);
  const int nz = std::max(1, count_unique(zs) - 1);
  return {nx, ny, nz};
}

// Extend a box mesh upward by `thickness`; preserve the original nx/ny but
// add `nz_add` new layers in the extension region only.
// Used by photo() (coarse is acceptable) — deposit() uses extend_mesh_exact().
Mesh extend_mesh(const Mesh& base, double thickness, int nz_add) {
  const BBox bb = base.bbox();
  const double base_lz = bb.hi.z - bb.lo.z;
  const int nc = static_cast<int>(base.cells.size());
  const int nx_eff = std::max(2, static_cast<int>(
      std::round(std::sqrt(std::sqrt(static_cast<double>(nc) / 6.0 *
                           (bb.hi.x - bb.lo.x) * (bb.hi.y - bb.lo.y) /
                           base_lz)))));
  const double cell_h = base_lz / std::max(1.0,
      std::round(base_lz * std::cbrt(nc / 6.0) / std::sqrt(
          (bb.hi.x - bb.lo.x) * (bb.hi.y - bb.lo.y))));
  const int nz_base_est = std::max(1, static_cast<int>(std::round(base_lz / cell_h)));
  const int nz_total = nz_base_est + nz_add;
  return make_box_mesh(bb.lo.x, bb.hi.x, bb.lo.y, bb.hi.y, bb.lo.z,
                       bb.lo.z + base_lz + thickness, nx_eff, nx_eff, nz_total);
}

// Extend a box mesh upward, preserving the original nx and ny exactly, and
// adding `nz_add` layers proportionally sized to existing z-cell height.
// Used by deposit() where mesh quality must be maintained.
Mesh extend_mesh_exact(const Mesh& base, double thickness, int nz_add) {
  const BBox bb = base.bbox();
  const auto [nx, ny, nz] = infer_box_dims(base);
  const double base_lz = bb.hi.z - bb.lo.z;
  const int nz_total = nz + nz_add;
  return make_box_mesh(bb.lo.x, bb.hi.x, bb.lo.y, bb.hi.y, bb.lo.z,
                       bb.lo.z + base_lz + thickness, nx, ny, nz_total);
}

const Dopant* dopant_or_throw(const std::string& species) {
  const Dopant* d = find_dopant(species);
  if (!d) throw std::runtime_error("unknown species '" + species + "' (B, P, As, Sb)");
  return d;
}

void need_mesh(const SimState& st) {
  if (!st.has_mesh) throw std::runtime_error("no mesh defined yet");
}

// "+1" damage model: an implant creates roughly one excess self-interstitial
// per implanted ion, distributed like the added dopant profile. Accumulate that
// excess into the "I" field (cm^-3), which diffuse_ted() consumes. `before` is
// the dopant field prior to this implant (empty means all-zero baseline).
void seed_interstitials(SimState& st, const std::vector<double>& before,
                        const std::vector<double>& after) {
  auto& I = st.fields["I"];
  I.resize(st.mesh.cells.size(), 0.0);
  const double survival = ParamDB::instance().get("ted.frenkel_survival", 1.0);
  for (std::size_t i = 0; i < after.size(); ++i) {
    const double added = after[i] - (i < before.size() ? before[i] : 0.0);
    if (added > 0) I[i] += survival * added;
  }
}

// Seed the interstitial excess from the MC Kinchin-Pease damage field:
//   I += kFrenkelSurvival * damage, capped at kAmorphizationDensity.
void seed_interstitials_from_damage(SimState& st,
                                    const std::vector<double>& damage) {
  auto& I = st.fields["I"];
  I.resize(st.mesh.cells.size(), 0.0);
  for (std::size_t i = 0; i < damage.size() && i < I.size(); ++i) {
    I[i] += kFrenkelSurvival * damage[i];
    I[i] = std::min(I[i], kAmorphizationDensity);
  }
}

}  // namespace

std::vector<char> silicon_mask(const SimState& st) {
  std::vector<char> mask(st.mesh.cells.size(), 1);
  for (std::size_t i = 0; i < st.mesh.cells.size(); ++i) {
    auto it = st.region_material.find(st.mesh.cell_region[i]);
    mask[i] = (it == st.region_material.end()) || is_silicon(it->second);
  }
  return mask;
}

// Full MatId per cell (Si/oxide/nitride/poly/gas); untagged cells are Si.
std::vector<int> material_ids(const SimState& st) {
  std::vector<int> mat(st.mesh.cells.size(), kMatSi);
  for (std::size_t i = 0; i < st.mesh.cells.size(); ++i) {
    auto it = st.region_material.find(st.mesh.cell_region[i]);
    mat[i] = (it == st.region_material.end()) ? static_cast<int>(kMatSi)
                                              : static_cast<int>(material_id(it->second));
  }
  return mat;
}

int resolve_region(const SimState& st, const std::string& v) {
  char* end = nullptr;
  const long tag = std::strtol(v.c_str(), &end, 10);
  if (end != v.c_str() && *end == '\0') return static_cast<int>(tag);
  for (const auto& [t, n] : st.mesh.region_names)
    if (lower(n) == lower(v)) return t;
  throw std::runtime_error("unknown region '" + v + "'");
}

void mesh_box(SimState& st, double x0, double x1, double y0, double y1,
              double z0, double z1, int nx, int ny, int nz, std::ostream* log) {
  st.mesh = make_box_mesh(x0, x1, y0, y1, z0, z1, nx, ny, nz);
  st.has_mesh = true;
  st.fields.clear();
  st.bcs.clear();
  st.region_material.clear();
  st.layer_stack.clear();
  strip(st, nullptr);
  for (int tag : st.mesh.region_tags()) st.region_material[tag] = "silicon";
  if (log) {
    const BBox b = st.mesh.bbox();
    *log << "[mesh] " << st.mesh.cells.size() << " tets, " << st.mesh.nodes.size()
         << " nodes, extent " << fmt("%.3g", (b.hi.x - b.lo.x) * 1e4) << " x "
         << fmt("%.3g", (b.hi.y - b.lo.y) * 1e4) << " x "
         << fmt("%.3g", (b.hi.z - b.lo.z) * 1e4) << " um, min orthogonality "
         << fmt("%.3f", st.mesh.min_orthogonality()) << "\n";
  }
}

void mesh_gmsh(SimState& st, const std::string& file, double scale,
               std::ostream* log) {
  st.mesh = read_gmsh(file, scale, log);
  st.has_mesh = true;
  st.fields.clear();
  st.bcs.clear();
  st.region_material.clear();
  st.layer_stack.clear();
  strip(st, nullptr);
  for (int tag : st.mesh.region_tags()) st.region_material[tag] = "silicon";
}

void set_region(SimState& st, const std::string& material, int tag,
                std::ostream* log) {
  need_mesh(st);
  const std::string mat = lower(material);
  static const char* known[] = {"silicon", "si", "oxide", "nitride", "poly",
                                "polysilicon", "gas"};
  if (std::none_of(std::begin(known), std::end(known),
                   [&](const char* k) { return mat == k; }))
    throw std::runtime_error("unknown material '" + material + "'");
  if (tag < 0) {
    for (int t : st.mesh.region_tags()) st.region_material[t] = mat;
  } else {
    st.region_material[tag] = mat;
  }
  if (log) *log << "[region] -> " << mat << "\n";
}

void init(SimState& st, const std::string& species, double conc, int region,
          std::ostream* log) {
  need_mesh(st);
  const Dopant* d = dopant_or_throw(species);
  if (conc < 0) throw std::runtime_error("conc must be >= 0");
  auto& f = st.fields[d->symbol];
  f.resize(st.mesh.cells.size(), 0.0);
  const std::vector<char> mask = silicon_mask(st);
  std::size_t nset = 0;
  for (std::size_t i = 0; i < f.size(); ++i) {
    if (!mask[i]) continue;
    if (region >= 0 && st.mesh.cell_region[i] != region) continue;
    f[i] = conc;
    ++nset;
  }
  if (log)
    *log << "[init] " << d->symbol << " = " << fmt("%.3g", conc) << " cm^-3 in "
         << nset << " cells\n";
}

double implant_gauss(SimState& st, const std::string& species, double dose,
                     double energy_kev, double rp, double drp, double drl,
                     bool has_window, double x1, double x2, double y1, double y2,
                     bool seed_damage, const std::string& profile,
                     std::ostream* log) {
  need_mesh(st);
  if (profile != "gauss" && profile != "pearson")
    throw std::runtime_error("implant: profile must be gauss|pearson");
  ImplantParams p;
  p.dopant = dopant_or_throw(species);
  p.dose = dose;
  if (profile == "pearson") {
    if (!(energy_kev > 0))
      throw std::runtime_error("pearson profile requires energy=");
    double gamma = 0, beta = 3;
    if (!implant_moments(*p.dopant, energy_kev, p.rp, p.drp, gamma, beta))
      throw std::runtime_error("no range table for species; give rp/drp");
    p.profile = ImplantParams::Profile::pearson4;
    p.gamma = gamma;
    p.beta = beta;
  } else if (energy_kev > 0) {
    if (!implant_range(*p.dopant, energy_kev, p.rp, p.drp))
      throw std::runtime_error("no range table for species; give rp/drp");
  } else {
    p.rp = rp;
    p.drp = drp;
  }
  p.drl = drl;
  p.has_window = has_window;
  p.x1 = x1; p.x2 = x2; p.y1 = y1; p.y2 = y2;
  auto& f = st.fields[p.dopant->symbol];
  f.resize(st.mesh.cells.size(), 0.0);
  const std::vector<double> before = seed_damage ? f : std::vector<double>{};
  const double atoms = apply_implant(st.mesh, silicon_mask(st), p, f);
  if (seed_damage) seed_interstitials(st, before, f);
  bool used_pearson = false;
  if (p.profile == ImplantParams::Profile::pearson4) {
    // Mirror apply_implant's Type-IV validity check for accurate logging
    // (apply_implant itself falls back silently to Gaussian when invalid).
    const double A = 10.0 * p.beta - 12.0 * p.gamma * p.gamma - 18.0;
    if (A != 0.0 && p.beta > 1.0 + p.gamma * p.gamma) {
      const double b0 = -p.drp * p.drp * (4.0 * p.beta - 3.0 * p.gamma * p.gamma) / A;
      const double b1 = -p.gamma * p.drp * (p.beta + 3.0) / A;
      const double b2 = -(2.0 * p.beta - 3.0 * p.gamma * p.gamma - 6.0) / A;
      used_pearson = (b1 * b1 - 4.0 * b0 * b2) < 0.0;
    }
  }
  if (log)
    *log << "[implant] " << p.dopant->symbol << " gauss: dose="
         << fmt("%.3g", p.dose) << " cm^-2, Rp=" << fmt("%.4g", p.rp * 1e4)
         << " um, dRp=" << fmt("%.4g", p.drp * 1e4)
         << " um, profile=" << (used_pearson ? "pearson4" : "gauss") << "\n";
  return atoms;
}

McImplantStats implant_mc(SimState& st, const std::string& species, double dose,
                          double energy_kev, long long ions, double tilt_deg,
                          double rotation_deg, unsigned long long seed,
                          int threads, bool channeling, bool has_window,
                          double x1, double x2, double y1, double y2,
                          bool seed_damage, std::ostream* log) {
  need_mesh(st);
  const Dopant* dop = dopant_or_throw(species);
  auto& f = st.fields[dop->symbol];
  f.resize(st.mesh.cells.size(), 0.0);
  const std::vector<double> before = seed_damage ? f : std::vector<double>{};
  const bool use_damage = seed_damage && channeling;

  McImplantParams p;
  p.dopant = dop;
  p.dose = dose;
  p.energy_kev = energy_kev;
  p.tilt_deg = tilt_deg;
  p.rotation_deg = rotation_deg;
  p.ions = ions;
  p.threads = threads;
  p.seed = seed;
  p.channeling = channeling;
  p.has_window = has_window;
  p.x1 = x1; p.x2 = x2; p.y1 = y1; p.y2 = y2;

  // Physical resist stack present: transport through the full stack and
  // transfer the silicon profile back onto the working mesh.
  if (st.has_stack) {
    const int nc_s = static_cast<int>(st.stack.cells.size());
    std::vector<char> stack_si(nc_s, 0);
    for (int ci = 0; ci < nc_s; ++ci)
      stack_si[ci] = (st.stack_cell_mat[ci] == 0) ? 1 : 0;

    p.has_window = false;  // full-surface beam; masking is physical
    p.material_table = st.mat_table;
    p.cell_material = &st.stack_cell_mat;

    std::vector<double> stack_conc(nc_s, 0.0);
    std::vector<double> stack_dmg;
    const McImplantStats s = apply_mc_implant(
        st.stack, stack_si, p, stack_conc, use_damage ? &stack_dmg : nullptr);

    const std::vector<double> transferred =
        transfer_field_nearest(st.stack, stack_conc, st.mesh);
    for (std::size_t i = 0; i < f.size(); ++i) f[i] += transferred[i];
    if (use_damage) {
      const std::vector<double> dmg_transferred =
          transfer_field_nearest(st.stack, stack_dmg, st.mesh);
      seed_interstitials_from_damage(st, dmg_transferred);
    } else if (seed_damage) {
      seed_interstitials(st, before, f);
    }

    if (log) {
      const double n = static_cast<double>(p.ions);
      *log << "[implant] " << dop->symbol << " MC (physical resist): E="
           << fmt("%.4g", p.energy_kev) << " keV, deposited_in_Si="
           << fmt("%.1f", 100.0 * s.deposited / n) << "%, stopped_in_resist="
           << fmt("%.1f", 100.0 * s.in_mask / n) << "%\n";
      if (seed_damage) {
        if (use_damage) {
          double peak = 0.0;
          for (double v : st.fields["I"]) peak = std::max(peak, v);
          *log << "[implant] damage seed: peak I=" << fmt("%.3g", peak)
               << " cm^-3 (KP damage x " << fmt("%.3g", kFrenkelSurvival)
               << ")\n";
        } else {
          *log << "[implant] damage: '+1' model (channeling off)\n";
        }
      }
    }
    return s;
  }

  std::vector<double> dmg;
  const McImplantStats s = apply_mc_implant(st.mesh, silicon_mask(st), p, f,
                                            use_damage ? &dmg : nullptr);
  if (use_damage) {
    seed_interstitials_from_damage(st, dmg);
  } else if (seed_damage) {
    seed_interstitials(st, before, f);
  }
  if (log) {
    const double n = static_cast<double>(p.ions);
    *log << "[implant] " << dop->symbol << " MC: E=" << fmt("%.4g", p.energy_kev)
         << " keV, deposited " << fmt("%.1f", 100.0 * s.deposited / n)
         << "% (Rp=" << fmt("%.4g", s.rp * 1e4) << " um)\n";
    if (seed_damage) {
      if (use_damage) {
        double peak = 0.0;
        for (double v : st.fields["I"]) peak = std::max(peak, v);
        *log << "[implant] damage seed: peak I=" << fmt("%.3g", peak)
             << " cm^-3 (KP damage x " << fmt("%.3g", kFrenkelSurvival)
             << ")\n";
      } else {
        *log << "[implant] damage: '+1' model (channeling off)\n";
      }
    }
  }
  return s;
}

void photo(SimState& st, double thickness, int nz_add, std::ostream* log) {
  need_mesh(st);
  if (st.has_stack) throw std::runtime_error("resist already present; strip first");
  if (thickness <= 0) throw std::runtime_error("resist thickness must be > 0");

  const BBox bb = st.mesh.bbox();
  const double z_si_top = bb.hi.z;
  st.stack = extend_mesh(st.mesh, thickness, nz_add);
  const int nc_s = static_cast<int>(st.stack.cells.size());

  st.stack_cell_mat.assign(nc_s, 0);
  for (int ci = 0; ci < nc_s; ++ci)
    if (st.stack.cell_cent[ci].z > z_si_top) st.stack_cell_mat[ci] = 1;  // resist

  st.mat_table = {target_silicon(), target_photoresist(), target_vacuum()};
  st.stack_resist_z0 = z_si_top;
  st.has_stack = true;

  if (log) {
    long rc = 0;
    for (int m : st.stack_cell_mat) if (m == 1) ++rc;
    *log << "[photo] resist=" << fmt("%.4g", thickness * 1e4) << " um, stack "
         << nc_s << " tets, " << rc << " resist cells\n";
  }
}

// Winding-number point-in-polygon test for a closed polygon given as (x,y) pairs.
// Returns true when (px,py) is strictly inside or on the boundary.
static bool point_in_polygon(double px, double py,
                              const std::vector<std::pair<double,double>>& poly) {
  const int n = static_cast<int>(poly.size());
  if (n < 3) return false;
  int winding = 0;
  for (int i = 0; i < n; ++i) {
    const auto [ax, ay] = poly[i];
    const auto [bx, by] = poly[(i + 1) % n];
    if (ay <= py) {
      if (by > py) {
        // upward crossing
        const double cross = (bx - ax) * (py - ay) - (by - ay) * (px - ax);
        if (cross > 0) ++winding;
      }
    } else {
      if (by <= py) {
        // downward crossing
        const double cross = (bx - ax) * (py - ay) - (by - ay) * (px - ax);
        if (cross < 0) --winding;
      }
    }
  }
  return winding != 0;
}

void mask(SimState& st, double x1, double x2, double y1, double y2,
          std::ostream* log) {
  need_mesh(st);
  if (!st.has_stack) throw std::runtime_error("no resist present; photo first");
  if (!(x2 > x1)) throw std::runtime_error("mask requires x2 > x1");
  if (!(y2 > y1)) throw std::runtime_error("mask requires y2 > y1");
  // Delegate to mask_polygon with a rectangle.
  mask_polygon(st, {{x1, y1}, {x2, y1}, {x2, y2}, {x1, y2}}, log);
}

void mask_polygon(SimState& st,
                  const std::vector<std::pair<double,double>>& poly,
                  std::ostream* log) {
  need_mesh(st);
  if (!st.has_stack) throw std::runtime_error("no resist present; photo first");
  if (poly.size() < 3) throw std::runtime_error("mask_polygon requires >= 3 vertices");
  int opened = 0;
  const int nc_s = static_cast<int>(st.stack.cells.size());
  for (int ci = 0; ci < nc_s; ++ci) {
    if (st.stack_cell_mat[ci] != 1) continue;
    const Vec3& c = st.stack.cell_cent[ci];
    if (point_in_polygon(c.x, c.y, poly)) {
      st.stack_cell_mat[ci] = 2;  // vacuum / developed opening
      ++opened;
    }
  }
  if (log)
    *log << "[mask_polygon] " << poly.size() << " vertices -> opened "
         << opened << " cells\n";
}

void strip(SimState& st, std::ostream* log) {
  if (!st.has_stack) {
    if (log) *log << "[strip] no resist present\n";
    return;
  }
  st.has_stack = false;
  st.stack = Mesh{};
  st.stack_cell_mat.clear();
  st.mat_table.clear();
  st.stack_resist_z0 = 0;
  if (log) *log << "[strip] photoresist removed\n";
}

void deposit(SimState& st, const std::string& material,
             double thickness, int nz_add,
             const std::vector<std::pair<double,double>>& poly,
             std::ostream* log) {
  need_mesh(st);
  if (thickness <= 0) throw std::runtime_error("deposit: thickness must be > 0");

  // Validate material name.
  static const char* known[] = {"oxide", "sio2", "nitride", "si3n4", "poly",
                                 "polysilicon", "silicon", "si"};
  const std::string mat = lower(material);
  if (std::none_of(std::begin(known), std::end(known),
                   [&](const char* k) { return mat == k; }))
    throw std::runtime_error("deposit: unknown material '" + material + "'");

  // Extend the mesh upward preserving original nx/ny resolution.
  Mesh ext = extend_mesh_exact(st.mesh, thickness, nz_add);
  const BBox bb = st.mesh.bbox();
  const double z_top = bb.hi.z;

  // Tag newly added cells as `material` when they are (a) above the old surface
  // and (b) inside the polygon (or unconditionally if poly is empty).
  const bool use_poly = poly.size() >= 3;
  // Bug fix (P1-7): only default *unknown* tags to "silicon" — a blanket
  // "silicon" write here used to clobber the region_material of prior
  // deposited films (e.g. a first deposit("oxide", ...) layer) whose tag
  // survives into `ext.region_tags()` unchanged.
  for (int tag : ext.region_tags())
    if (st.region_material.find(tag) == st.region_material.end())
      st.region_material[tag] = "silicon";

  // We store the deposit material in the region_material map keyed by a new
  // synthetic region tag that doesn't conflict with existing ones.
  const auto existing_tags = st.mesh.region_tags();
  const int max_tag = existing_tags.empty() ? 0
      : *std::max_element(existing_tags.begin(), existing_tags.end());
  const int dep_tag = max_tag + 1000;  // synthetic tag for deposited film

  const int nc_ext = static_cast<int>(ext.cells.size());
  ext.cell_region.resize(nc_ext, 0);
  // Keep existing cell_region for the base cells; set deposited cells.
  // ext is a fresh mesh, so re-populate cell_region from centroids.
  const int nc_old_r = static_cast<int>(st.mesh.cells.size());
  for (int ci = 0; ci < nc_ext; ++ci) {
    const Vec3& c = ext.cell_cent[ci];
    if (c.z > z_top) {
      // New deposited cell: apply polygon filter.
      if (!use_poly || point_in_polygon(c.x, c.y, poly))
        ext.cell_region[ci] = dep_tag;
    } else {
      // Bug fix (P1-7): base cells (z <= z_top) must inherit the region tag
      // of their nearest-centroid old cell, not the default-initialized 0 —
      // otherwise a prior film's material tag is lost after the next
      // deposit's mesh rebuild.
      double best = 1e300;
      int best_j = 0;
      for (int j = 0; j < nc_old_r; ++j) {
        const Vec3& oc = st.mesh.cell_cent[j];
        const double dx = c.x - oc.x, dy = c.y - oc.y, dz = c.z - oc.z;
        const double d2 = dx*dx + dy*dy + dz*dz;
        if (d2 < best) { best = d2; best_j = j; }
      }
      ext.cell_region[ci] = st.mesh.cell_region[best_j];
    }
  }
  st.region_material[dep_tag] = mat;

  // Transfer all existing fields to the extended mesh. Cells above the old
  // outer surface belong to the freshly deposited film and start clean (zero
  // concentration) — nearest-centroid transfer must not copy the substrate's
  // near-surface tail into them.
  for (auto& [sym, conc] : st.fields) {
    std::vector<double> new_conc(nc_ext, 0.0);
    for (int ci = 0; ci < nc_ext; ++ci) {
      const Vec3& cc = ext.cell_cent[ci];
      if (cc.z > z_top) continue;  // new film cell: stays 0
      double best = 1e300;
      int best_j = 0;
      const int nc_old = static_cast<int>(st.mesh.cells.size());
      for (int j = 0; j < nc_old; ++j) {
        const Vec3& oc = st.mesh.cell_cent[j];
        const double dx = cc.x - oc.x, dy = cc.y - oc.y, dz = cc.z - oc.z;
        const double d2 = dx*dx + dy*dy + dz*dz;
        if (d2 < best) { best = d2; best_j = j; }
      }
      new_conc[ci] = conc[best_j];
    }
    conc = std::move(new_conc);
  }

  st.mesh = std::move(ext);
  st.layer_stack.insert(st.layer_stack.begin(), {dep_tag, mat});

  if (log)
    *log << "[deposit] " << material << " thickness="
         << fmt("%.4g", thickness * 1e4) << " um"
         << (use_poly ? " (polygon masked)" : " (blanket)")
         << ", mesh now " << st.mesh.cells.size() << " tets\n";
}

namespace {
// True if `cell_mat` should be treated as etched under selectivity filter
// `filter` (already lower-cased is not required from caller). Empty filter
// matches every non-gas material. Non-empty filter matches case-insensitive,
// with the silicon "si"/"silicon" alias (same normalization as is_silicon).
bool material_matches(const std::string& cell_mat, const std::string& filter) {
  const std::string cm = lower(cell_mat);
  if (filter.empty()) return cm != "gas";
  const std::string fm = lower(filter);
  if (is_silicon(fm)) return is_silicon(cm);
  return cm == fm;
}
}  // namespace

void etch(SimState& st, double depth,
          const std::vector<std::pair<double,double>>& poly,
          const std::string& material,
          std::ostream* log) {
  need_mesh(st);
  if (depth <= 0) throw std::runtime_error("etch: depth must be > 0");
  if (st.has_stack) throw std::runtime_error("etch: strip resist first");

  const BBox bb = st.mesh.bbox();
  const double z_top = bb.hi.z;
  const double z_cut = z_top - depth;
  const bool use_poly = poly.size() >= 3;

  auto material_of = [&](int tag) -> std::string {
    auto it = st.region_material.find(tag);
    return it == st.region_material.end() ? "silicon" : it->second;
  };

  const int nc = static_cast<int>(st.mesh.cells.size());
  if (st.mesh.cell_region.size() != static_cast<std::size_t>(nc))
    st.mesh.cell_region.resize(nc, 0);

  if (!use_poly) {
    // ── Blanket etch: TRUE cell removal (node compaction + identity field
    // transfer). See process.hpp for why polygon etch keeps the old
    // gas-retag behavior instead.
    std::vector<char> keep(nc, 1);
    int removed = 0;
    for (int ci = 0; ci < nc; ++ci) {
      const Vec3& c = st.mesh.cell_cent[ci];
      if (c.z > z_cut &&
          material_matches(material_of(st.mesh.cell_region[ci]), material)) {
        keep[ci] = 0;
        ++removed;
      }
    }

    // Map old cell index -> new cell index (kept cells only).
    std::vector<int> new_cell_of(nc, -1);
    std::vector<int> kept_old;
    kept_old.reserve(nc - removed);
    for (int ci = 0; ci < nc; ++ci) {
      if (!keep[ci]) continue;
      new_cell_of[ci] = static_cast<int>(kept_old.size());
      kept_old.push_back(ci);
    }

    // Node compaction: collect nodes referenced by kept cells only.
    const int nn_old = static_cast<int>(st.mesh.nodes.size());
    std::vector<int> new_node_of(nn_old, -1);
    Mesh nm;
    nm.nodes.reserve(nn_old);
    for (int oci : kept_old) {
      for (int k = 0; k < 4; ++k) {
        const int nid = st.mesh.cells[oci][k];
        if (new_node_of[nid] < 0) {
          new_node_of[nid] = static_cast<int>(nm.nodes.size());
          nm.nodes.push_back(st.mesh.nodes[nid]);
        }
      }
    }
    nm.cells.resize(kept_old.size());
    nm.cell_region.resize(kept_old.size());
    for (std::size_t k = 0; k < kept_old.size(); ++k) {
      const int oci = kept_old[k];
      for (int v = 0; v < 4; ++v)
        nm.cells[k][v] = new_node_of[st.mesh.cells[oci][v]];
      nm.cell_region[k] = st.mesh.cell_region[oci];
    }
    nm.patch_names = st.mesh.patch_names;
    nm.region_names = st.mesh.region_names;
    nm.finalize();

    // Identity field transfer: kept cells copy their old value exactly, no
    // interpolation — mass in surviving cells is preserved bit-for-bit.
    for (auto& [sym, conc] : st.fields) {
      std::vector<double> new_conc(kept_old.size());
      for (std::size_t k = 0; k < kept_old.size(); ++k)
        new_conc[k] = conc[kept_old[k]];
      conc = std::move(new_conc);
    }

    st.mesh = std::move(nm);

    // Drop layer_stack entries for any region tag that has no surviving
    // cells (region_material entries themselves are left intact).
    auto tag_has_cells = [&](int tag) {
      for (int r : st.mesh.cell_region)
        if (r == tag) return true;
      return false;
    };
    for (auto it = st.layer_stack.begin(); it != st.layer_stack.end();) {
      if (!tag_has_cells(it->first)) it = st.layer_stack.erase(it);
      else ++it;
    }

    if (log)
      *log << "[etch] depth=" << fmt("%.4g", depth * 1e4) << " um (blanket"
           << (material.empty() ? "" : ", material=" + material)
           << "), removed " << removed << " cells, mesh now "
           << st.mesh.cells.size() << " tets\n";
    return;
  }

  // ── Polygon etch: gas re-tag (mesh topology unchanged).
  int gas_tag = -1;
  for (const auto& [t, m] : st.region_material)
    if (m == "gas") { gas_tag = t; break; }
  if (gas_tag < 0) {
    const auto etags = st.mesh.region_tags();
    const int max_tag = etags.empty() ? 0
        : *std::max_element(etags.begin(), etags.end());
    gas_tag = max_tag + 2000;
    st.region_material[gas_tag] = "gas";
  }

  int etched = 0;
  for (int ci = 0; ci < nc; ++ci) {
    const Vec3& c = st.mesh.cell_cent[ci];
    if (c.z < z_cut) continue;  // below etch depth
    if (!point_in_polygon(c.x, c.y, poly)) continue;
    if (!material_matches(material_of(st.mesh.cell_region[ci]), material)) continue;
    st.mesh.cell_region[ci] = gas_tag;
    ++etched;
  }

  // Zero out concentrations in etched cells (material removed).
  for (auto& [sym, conc] : st.fields)
    for (int ci = 0; ci < nc; ++ci)
      if (st.mesh.cell_region[ci] == gas_tag) conc[ci] = 0.0;

  if (log)
    *log << "[etch] depth=" << fmt("%.4g", depth * 1e4) << " um (polygon masked"
         << (material.empty() ? "" : ", material=" + material)
         << "), etched " << etched << " cells\n";
}

namespace {

// Advance `phi` by `total_time` under nodal velocity `F`, alternating
// short advect chunks (~5 CFL substeps each, per P2-5 step (c)) with
// `levelset_reinit` calls to keep it close to a true signed distance field
// (the advection itself only enforces the HJ transport equation, not the
// |∇φ|=1 eikonal constraint -- without periodic reinit the field would
// drift and the node-gradient-based velocity/normal computations would
// degrade over many substeps).
void advect_with_reinit(const Mesh& m, const MeshTopology& topo,
                        std::vector<double>& phi, const std::vector<double>& F,
                        double total_time, std::ostream* /*log*/) {
  double max_F = 0.0;
  for (double f : F) max_F = std::max(max_F, std::fabs(f));
  if (max_F <= 0.0 || total_time <= 0.0) return;

  const int nn = static_cast<int>(m.nodes.size());
  double h_min = std::numeric_limits<double>::max();
  for (int i = 0; i < nn; ++i)
    for (int j : topo.node_adj[i]) {
      const Vec3 dv = m.nodes[j] - m.nodes[i];
      const double len = std::sqrt(dot(dv, dv));
      if (len > 1e-30) h_min = std::min(h_min, len);
    }
  if (h_min == std::numeric_limits<double>::max()) return;

  const double dt_sub = 0.5 * h_min / max_F;
  // Reinit cadence and iteration count were tuned against two conflicting
  // standalone-probe measurements rather than following the spec's literal
  // "advect 5 substeps, then levelset_reinit" suggestion verbatim:
  //  - Lumping many pseudo-time iterations into one levelset_reinit() call
  //    made after a *chunk* of several CFL substeps (i.e. reinit(5) called
  //    once every 5 substeps, as literally suggested) is unstable at the
  //    mesh resolutions these acceptance tests run at: levelset_reinit's
  //    explicit Sussman-PDE update let sign flips appear tens of node-hops
  //    away from the true front in the isotropic mask-undercut scenario
  //    (measured: max_iters 2/3/4/5 in that one lumped call -> 48/64/96/112
  //    spuriously "etched" cells far past the physically-reachable front).
  //    Reinitializing every *single* CFL substep instead (finer cadence,
  //    same total pseudo-time budget) is far more stable.
  //  - But the vertical/anisotropic etch-depth test (acceptance criterion
  //    1, sensitive to phi's *magnitude* via the avg-of-4-node-phi > 0
  //    retag threshold, not just its sign) needs enough redistancing to
  //    converge |grad phi| -> 1 well enough to land within one cell height
  //    of the geometric etch(depth=) reference: max_iters=1 or 2 per
  //    substep under-corrects and the retag threshold is never crossed
  //    (measured depth stuck at 0 cells removed); max_iters=4 per substep
  //    hits the reference exactly, and (measured against the same
  //    mask-undercut scenario used above) is still far short of the
  //    instability onset seen with the lumped cadence -- so 4 substep-local
  //    iterations is the value used here, chosen as the best point
  //    satisfying both measurements simultaneously.
  double remaining = total_time;
  while (remaining > 1e-30) {
    const double step = std::min(dt_sub, remaining);
    levelset_advect(m, topo, phi, F, step);
    remaining -= step;
    levelset_reinit(m, phi, 4);
  }
}

// Extend a per-node etch-rate "seed" (valid only at true solid/gas interface
// nodes) to the whole mesh by nearest-interface propagation (Dijkstra over
// mesh edges, mirroring levelset_init's distance BFS). This matters whenever
// a slow/no-etch material sits directly against a fast-etching one *without*
// itself bordering gas (e.g. a thin nitride mask sitting on exposed
// silicon): naively taking "the fastest rate among any incident solid cell"
// at every node (as the spec's plain-language description suggests) also
// paints the *buried* silicon/nitride contact plane with the silicon rate,
// even though that plane never touches gas -- found via a standalone probe
// reproducing the isotropic-undercut acceptance test, where it caused the
// entire masked region to erode near-instantly instead of a bounded lateral
// undercut. Propagating each node's rate from its nearest true solid/gas
// interface point (rather than from any nearby solid material) fixes this:
// points buried under the mask correctly inherit the mask's own (zero) rate
// from the interface directly above them, and only pick up the substrate's
// rate once they are close enough to the actual opening for that to be the
// nearer interface -- which is exactly the physical lateral-undercut length
// scale. `inside_mask`: per-cell, matches levelset_init's convention.
// `node_seed_rate`: per-node, the naive "max incident inside-cell rate"
// (used only at nodes that do touch gas -- the true interface -- to seed
// the propagation). Returns the per-node extended rate field, and (via
// `true_interface_out` if non-null) which nodes are true solid/gas interface
// nodes (touch both an inside and an outside/gas cell) -- callers that also
// want an anisotropic (surface-normal) factor should restrict it to these
// nodes, for the same buried-interface reason.
std::vector<double> extend_velocity_from_interface(
    const Mesh& m, const MeshTopology& topo,
    const std::vector<char>& inside_mask,
    const std::vector<double>& node_seed_rate,
    std::vector<char>* true_interface_out = nullptr) {
  const int nn = static_cast<int>(m.nodes.size());
  const int nc = static_cast<int>(m.cells.size());
  std::vector<char> node_in(nn, 0), node_out(nn, 0);
  for (int ci = 0; ci < nc; ++ci) {
    const bool in = inside_mask[ci] != 0;
    for (int v : m.cells[ci]) { if (in) node_in[v] = 1; else node_out[v] = 1; }
  }

  std::vector<double> Fext(nn, 0.0);
  std::vector<double> dist(nn, std::numeric_limits<double>::max());
  using P = std::pair<double, int>;
  std::priority_queue<P, std::vector<P>, std::greater<P>> pq;
  std::vector<char> true_interface(nn, 0);
  for (int i = 0; i < nn; ++i) {
    if (node_in[i] && node_out[i]) {
      true_interface[i] = 1;
      dist[i] = 0.0;
      Fext[i] = node_seed_rate[i];
      pq.push({0.0, i});
    }
  }
  while (!pq.empty()) {
    const auto [d, u] = pq.top(); pq.pop();
    if (d > dist[u]) continue;
    for (int v : topo.node_adj[u]) {
      const Vec3 dv = m.nodes[v] - m.nodes[u];
      const double len = std::sqrt(dot(dv, dv));
      const double nd = dist[u] + len;
      if (nd < dist[v]) {
        dist[v] = nd;
        Fext[v] = Fext[u];
        pq.push({nd, v});
      }
    }
  }
  if (true_interface_out) *true_interface_out = std::move(true_interface);
  return Fext;
}

// Common headroom-extension step used by both etch_rate and
// deposit_conformal: grow the mesh upward by `headroom_cm` (tagged `tag`,
// via extend_mesh_exact's regular-box machinery), inherit region tags for
// the pre-existing cells by nearest centroid, and carry all fields across
// (identity for old cells, zero for the new headroom). Returns the new
// mesh; `z_top0` is the old top-of-mesh z (cells above it are new).
Mesh add_headroom(SimState& st, double headroom_cm, int new_tag, double& z_top0_out) {
  const auto [nx0, ny0, nz0] = infer_box_dims(st.mesh);
  const BBox bb0 = st.mesh.bbox();
  const double base_lz = bb0.hi.z - bb0.lo.z;
  const double h = base_lz / std::max(1, nz0);
  const int nz_add = std::max(1, static_cast<int>(std::ceil(headroom_cm / h)));
  // Snap the added thickness to an exact multiple of the existing cell
  // height h: extend_mesh_exact rebuilds a *uniform* grid over
  // [z0, z0+base_lz+thickness] with nz+nz_add layers, so unless
  // thickness == nz_add*h exactly, the new layer height drifts away from h
  // and the old top-of-mesh z (z_top0) no longer lands on a grid line --
  // producing a jagged (off-by-a-fraction-of-a-cell) solid/gas boundary
  // instead of the clean planar interface the level-set init/advection
  // (and compute_node_normals) assume. Found via a standalone probe: without
  // this snap, no mesh node existed exactly at z_top0 and etch_rate silently
  // etched nothing.
  const double headroom_snapped = nz_add * h;
  const double z_top0 = bb0.hi.z;
  z_top0_out = z_top0;

  Mesh ext = extend_mesh_exact(st.mesh, headroom_snapped, nz_add);
  const int nc_ext = static_cast<int>(ext.cells.size());
  ext.cell_region.resize(nc_ext, 0);
  const int nc_old = static_cast<int>(st.mesh.cells.size());
  for (int ci = 0; ci < nc_ext; ++ci) {
    const Vec3& c = ext.cell_cent[ci];
    if (c.z > z_top0) {
      ext.cell_region[ci] = new_tag;
      continue;
    }
    double best = 1e300;
    int best_j = 0;
    for (int j = 0; j < nc_old; ++j) {
      const Vec3& oc = st.mesh.cell_cent[j];
      const double dx = c.x - oc.x, dy = c.y - oc.y, dz = c.z - oc.z;
      const double d2 = dx * dx + dy * dy + dz * dz;
      if (d2 < best) { best = d2; best_j = j; }
    }
    ext.cell_region[ci] = st.mesh.cell_region[best_j];
  }

  for (auto& [sym, conc] : st.fields) {
    std::vector<double> new_conc(nc_ext, 0.0);
    for (int ci = 0; ci < nc_ext; ++ci) {
      const Vec3& cc = ext.cell_cent[ci];
      if (cc.z > z_top0) continue;  // new cell: stays 0
      double best = 1e300;
      int best_j = 0;
      for (int j = 0; j < nc_old; ++j) {
        const Vec3& oc = st.mesh.cell_cent[j];
        const double dx = cc.x - oc.x, dy = cc.y - oc.y, dz = cc.z - oc.z;
        const double d2 = dx * dx + dy * dy + dz * dz;
        if (d2 < best) { best = d2; best_j = j; }
      }
      new_conc[ci] = conc[best_j];
    }
    conc = std::move(new_conc);
  }
  return ext;
}

int gas_tag_of(SimState& st, int mint_offset) {
  for (const auto& [t, m] : st.region_material)
    if (m == "gas") return t;
  const auto etags = st.mesh.region_tags();
  const int max_tag = etags.empty() ? 0 : *std::max_element(etags.begin(), etags.end());
  const int gas_tag = max_tag + mint_offset;
  st.region_material[gas_tag] = "gas";
  return gas_tag;
}

}  // namespace

void etch_rate(SimState& st, const std::map<std::string, double>& rates,
              double time_s, bool isotropic,
              const std::vector<std::pair<double,double>>& poly,
              std::ostream* log) {
  need_mesh(st);
  if (time_s <= 0) throw std::runtime_error("etch_rate: time_s must be > 0");
  if (st.has_stack) throw std::runtime_error("etch_rate: strip resist first");
  const bool use_poly = poly.size() >= 3;

  auto material_of = [&](int tag) -> std::string {
    auto it = st.region_material.find(tag);
    return it == st.region_material.end() ? "silicon" : it->second;
  };
  auto rate_of = [&](const std::string& mat) -> double {
    const std::string m = lower(mat);
    if (m == "gas") return 0.0;
    for (const auto& [name, r] : rates) {
      const std::string rn = lower(name);
      if (rn == m || (is_silicon(rn) && is_silicon(m))) return r;
    }
    return 0.0;
  };

  double max_rate = 0.0;
  for (const auto& [name, r] : rates) max_rate = std::max(max_rate, r);

  const BBox bb_pre = st.mesh.bbox();
  const auto [nx_pre, ny_pre, nz_pre] = infer_box_dims(st.mesh);
  const double h_pre = (bb_pre.hi.z - bb_pre.lo.z) / std::max(1, nz_pre);
  // Headroom: at least one cell layer even if nothing will actually etch, so
  // there is always an explicit "gas" side for the level set to advect into
  // -- a fresh solid box has no exposed gas cells yet (spec: "gas cells (or
  // the top boundary)"); we always materialize the latter as the former.
  const double headroom = std::max(max_rate * time_s * 1.25, h_pre);

  const int gas_tag = gas_tag_of(st, 3000);
  double z_top0 = 0.0;
  Mesh ext = add_headroom(st, headroom, gas_tag, z_top0);

  const int nc = static_cast<int>(ext.cells.size());
  std::vector<char> inside_mask(nc, 0);
  for (int ci = 0; ci < nc; ++ci)
    inside_mask[ci] = (material_of(ext.cell_region[ci]) != "gas") ? 1 : 0;

  MeshTopology topo;
  topo.build(ext);
  std::vector<double> phi = levelset_init(ext, inside_mask);

  const int nn = static_cast<int>(ext.nodes.size());
  // Seed rate: naive "max incident inside-cell rate" per node, but only
  // *used* at nodes that turn out to be true solid/gas interface nodes (see
  // extend_velocity_from_interface) -- this is what avoids painting a
  // buried multi-material contact (e.g. a nitride mask sitting directly on
  // exposed silicon) with the substrate's rate.
  std::vector<double> node_seed_rate(nn, 0.0);
  for (int ci = 0; ci < nc; ++ci) {
    if (!inside_mask[ci]) continue;
    const Vec3& cc = ext.cell_cent[ci];
    if (use_poly && !point_in_polygon(cc.x, cc.y, poly)) continue;
    const double r = rate_of(material_of(ext.cell_region[ci]));
    if (r <= 0) continue;
    for (int v : ext.cells[ci]) node_seed_rate[v] = std::max(node_seed_rate[v], r);
  }
  std::vector<char> true_interface;
  std::vector<double> F =
      extend_velocity_from_interface(ext, topo, inside_mask, node_seed_rate, &true_interface);

  if (!isotropic) {
    // Anisotropic (vertical/RIE) etch: only remove material at *true*
    // solid/gas interface nodes (true_interface, not just any multi-region
    // interface -- see extend_velocity_from_interface) whose surface normal
    // (compute_node_normals, shared with ale_mover.cpp) has a *positive* z
    // component, i.e. an upward-facing boundary (n points from solid into
    // gas at the exposed top surface). Verified with a standalone probe
    // against a flat Si box: this sign etches the top surface downward as
    // expected; the opposite sign instead affected only the domain-floor
    // boundary nodes (whose outward normal points in -z).
    const auto normals = compute_node_normals(ext, topo);
    for (int i = 0; i < nn; ++i)
      F[i] = true_interface[i] ? F[i] * std::max(0.0, normals[i].z) : 0.0;
  }
  // levelset_advect solves phi_t + F|grad phi| = 0 with "inside" (solid) at
  // phi<0; a positive F GROWS the inside region (that's what
  // deposit_conformal below wants, with F=+1 uniformly). Etching does the
  // opposite -- it shrinks the solid -- so the (non-negative) etch-rate
  // magnitude built above must be negated before advecting: verified with a
  // standalone probe that without this negation phi moves further negative
  // at the exposed surface (the solid grows) instead of etching it away.
  for (int i = 0; i < nn; ++i) F[i] = -F[i];

  advect_with_reinit(ext, topo, phi, F, time_s, log);

  int etched = 0;
  for (int ci = 0; ci < nc; ++ci) {
    if (!inside_mask[ci]) continue;
    // Defense in depth: a material with rate 0 (e.g. a mask) must never be
    // consumed even if phi drifted positive there (extend_velocity_from_
    // interface already prevents this in practice by giving such cells'
    // nodes rate 0, but a masked polygon or floating-point roundoff could
    // still leave a stray positive avg; this guard makes the "0 rate is
    // never removed" guarantee absolute rather than merely typical).
    if (rate_of(material_of(ext.cell_region[ci])) <= 0) continue;
    double avg = 0.0;
    for (int v : ext.cells[ci]) avg += phi[v];
    avg /= 4.0;
    if (avg > 0.0) { ext.cell_region[ci] = gas_tag; ++etched; }
  }
  for (auto& [sym, conc] : st.fields)
    for (int ci = 0; ci < nc; ++ci)
      if (ext.cell_region[ci] == gas_tag) conc[ci] = 0.0;

  {
    auto tag_has_cells = [&](int tag) {
      for (int r : ext.cell_region) if (r == tag) return true;
      return false;
    };
    for (auto it = st.layer_stack.begin(); it != st.layer_stack.end();) {
      if (!tag_has_cells(it->first)) it = st.layer_stack.erase(it);
      else ++it;
    }
  }

  st.mesh = std::move(ext);

  if (log)
    *log << "[etch_rate] " << (isotropic ? "isotropic" : "vertical") << " t="
         << fmt("%.4g", time_s) << " s" << (use_poly ? " (polygon masked)" : "")
         << ", retagged " << etched << " cells to gas (headroom "
         << fmt("%.4g", headroom * 1e4) << " um). If the undercut/etch front"
         << " looks blocky, run proc::refine first for finer lateral"
         << " mesh resolution near the interface.\n";
}

void deposit_conformal(SimState& st, const std::string& material,
                       double thickness_cm, std::ostream* log) {
  need_mesh(st);
  if (thickness_cm <= 0)
    throw std::runtime_error("deposit_conformal: thickness must be > 0");
  static const char* known[] = {"oxide", "sio2", "nitride", "si3n4", "poly",
                                 "polysilicon", "silicon", "si"};
  const std::string mat = lower(material);
  if (std::none_of(std::begin(known), std::end(known),
                   [&](const char* k) { return mat == k; }))
    throw std::runtime_error("deposit_conformal: unknown material '" + material + "'");

  auto material_of = [&](int tag) -> std::string {
    auto it = st.region_material.find(tag);
    return it == st.region_material.end() ? "silicon" : it->second;
  };

  const int gas_tag = gas_tag_of(st, 3000);
  double z_top0 = 0.0;
  // Headroom bigger than the film so the level set has room to reach full
  // thickness above the highest point of the (possibly stepped) surface;
  // the conformal growth still only deposits `thickness_cm` normal-to-
  // surface everywhere, including down any existing sidewall.
  Mesh ext = add_headroom(st, thickness_cm * 1.25, gas_tag, z_top0);

  const int nc = static_cast<int>(ext.cells.size());
  std::vector<char> inside_mask(nc, 0);
  for (int ci = 0; ci < nc; ++ci)
    inside_mask[ci] = (material_of(ext.cell_region[ci]) != "gas") ? 1 : 0;

  MeshTopology topo;
  topo.build(ext);
  std::vector<double> phi = levelset_init(ext, inside_mask);

  const int nn = static_cast<int>(ext.nodes.size());
  // Isotropic uniform growth: F=1 everywhere, advect for time=thickness_cm
  // so the interface advances by exactly thickness_cm normal-to-surface
  // (this is what gives sidewall/step coverage "for free").
  const std::vector<double> F(nn, 1.0);
  advect_with_reinit(ext, topo, phi, F, thickness_cm, log);

  const auto etags_now = ext.region_tags();
  const int max_tag_now = etags_now.empty() ? 0
      : *std::max_element(etags_now.begin(), etags_now.end());
  const int dep_tag = max_tag_now + 4000;
  st.region_material[dep_tag] = mat;

  int deposited = 0;
  for (int ci = 0; ci < nc; ++ci) {
    if (inside_mask[ci]) continue;  // already solid before this call
    double avg = 0.0;
    for (int v : ext.cells[ci]) avg += phi[v];
    avg /= 4.0;
    if (avg < 0.0) { ext.cell_region[ci] = dep_tag; ++deposited; }
  }

  st.mesh = std::move(ext);
  st.layer_stack.insert(st.layer_stack.begin(), {dep_tag, mat});

  if (log)
    *log << "[deposit_conformal] " << material << " thickness="
         << fmt("%.4g", thickness_cm * 1e4) << " um, conformal on "
         << deposited << " cells (including sidewalls)\n";
}

namespace {
bool is_oxide(const std::string& mat) {
  const std::string m = lower(mat);
  return m == "oxide" || m == "sio2";
}
}  // namespace

double oxidize(SimState& st, double time_s, double temp_k, bool wet,
               std::ostream* log) {
  need_mesh(st);
  if (st.has_stack)
    throw std::runtime_error("oxidize: resist stack is present; strip first");
  if (time_s <= 0) throw std::runtime_error("oxidize: time_s must be > 0");

  // (b) Measure existing oxide thickness.
  const int nc0 = static_cast<int>(st.mesh.cells.size());
  auto material_of = [&](int tag) -> std::string {
    auto it = st.region_material.find(tag);
    return it == st.region_material.end() ? "silicon" : it->second;
  };

  double z_si_top = -1e300;
  for (int ci = 0; ci < nc0; ++ci) {
    if (!is_silicon(material_of(st.mesh.cell_region[ci]))) continue;
    const auto& cell = st.mesh.cells[ci];
    for (int k = 0; k < 4; ++k)
      z_si_top = std::max(z_si_top, st.mesh.nodes[cell[k]].z);
  }
  const BBox bb0 = st.mesh.bbox();
  const double z_top = bb0.hi.z;
  if (z_si_top < -1e299) z_si_top = z_top;  // no silicon at all; degenerate

  for (int ci = 0; ci < nc0; ++ci) {
    if (st.mesh.cell_cent[ci].z <= z_si_top) continue;
    const std::string mat = material_of(st.mesh.cell_region[ci]);
    if (!is_oxide(mat))
      throw std::runtime_error(
          "oxidize: top surface is not bare Si or SiO2");
  }
  double x0_cm = z_top - z_si_top;
  if (x0_cm < 1e-9) x0_cm = 0.0;

  // P2-3 OED: `oed.theta` (ParamDB, default 0.01) is the fraction of the
  // Si atomic flux consumed at the growing interface that is assumed to
  // inject as excess self-interstitials into the Si just below the
  // interface. theta == 0 disables OED entirely and reproduces the P1-6
  // single-shot (geometry-only) behavior exactly -- kept as a separate
  // code path below so that behavior can never regress.
  const double theta = ParamDB::instance().get("oed.theta", 0.01);

  if (theta <= 0.0) {
    // ---- Legacy P1-6 path: one Deal-Grove step for the whole time_s. ----
    const double x_new_um = deal_grove_step(x0_cm * 1e4, time_s / 60.0,
                                            temp_k - 273.15, wet);
    const double dx_ox = x_new_um * 1e-4 - x0_cm;  // grown oxide, cm

    if (dx_ox <= 0) {
      if (log)
        *log << "[oxidize] " << (wet ? "wet " : "dry ") << fmt("%.6g", temp_k)
             << " K " << fmt("%.6g", time_s) << " s: no growth, tox="
             << fmt("%.4g", x_new_um) << " um\n";
      st.last_temp = temp_k;
      return x_new_um * 1e-4;
    }

    // (d) Volume bookkeeping and shape realization.
    const double dx_si = 0.44 * dx_ox;
    const double rise = 0.56 * dx_ox;
    const auto [nx0, ny0, nz0] = infer_box_dims(st.mesh);
    const double base_lz = bb0.hi.z - bb0.lo.z;
    const double h = base_lz / std::max(1, nz0);
    const int nz_add = std::max(1, static_cast<int>(std::round(rise / h)));
    Mesh ext = extend_mesh_exact(st.mesh, rise, nz_add);

    // (e) Retag: new Si/SiO2 interface height.
    const double z_if = z_si_top - dx_si;
    int ox_tag = -1;
    for (const auto& [t, m] : st.region_material)
      if (t >= 1000 && is_oxide(m)) { ox_tag = t; break; }
    if (ox_tag < 0) {
      const auto etags = st.mesh.region_tags();
      const int max_tag = etags.empty() ? 0
          : *std::max_element(etags.begin(), etags.end());
      ox_tag = max_tag + 3000;
    }
    st.region_material[ox_tag] = "oxide";

    const int nc_ext = static_cast<int>(ext.cells.size());
    ext.cell_region.resize(nc_ext, 0);
    for (int ci = 0; ci < nc_ext; ++ci) {
      const Vec3& c = ext.cell_cent[ci];
      if (c.z > z_if) {
        ext.cell_region[ci] = ox_tag;
      } else {
        // Nearest-centroid lookup of the old mesh's region tag.
        double best = 1e300;
        int best_j = 0;
        const int nc_old = static_cast<int>(st.mesh.cells.size());
        for (int j = 0; j < nc_old; ++j) {
          const Vec3& oc = st.mesh.cell_cent[j];
          const double dx = c.x - oc.x, dy = c.y - oc.y, dz = c.z - oc.z;
          const double d2 = dx*dx + dy*dy + dz*dz;
          if (d2 < best) { best = d2; best_j = j; }
        }
        ext.cell_region[ci] = st.mesh.cell_region[best_j];
      }
    }

    // (f) Field transfer: nearest-centroid, then zero dopants above the old
    // outer surface and keep the Si->SiO2 converted band frozen.
    for (auto& [sym, conc] : st.fields) {
      std::vector<double> new_conc(nc_ext, 0.0);
      for (int ci = 0; ci < nc_ext; ++ci) {
        const Vec3& cc = ext.cell_cent[ci];
        double best = 1e300;
        int best_j = 0;
        const int nc_old = static_cast<int>(st.mesh.cells.size());
        for (int j = 0; j < nc_old; ++j) {
          const Vec3& oc = st.mesh.cell_cent[j];
          const double dx = cc.x - oc.x, dy = cc.y - oc.y, dz = cc.z - oc.z;
          const double d2 = dx*dx + dy*dy + dz*dz;
          if (d2 < best) { best = d2; best_j = j; }
        }
        double v = conc[best_j];
        if (cc.z > z_si_top + x0_cm) v = 0.0;  // new cell above the old top
        new_conc[ci] = v;
      }
      conc = std::move(new_conc);
    }

    st.mesh = std::move(ext);

    // (h) Quality repair (extend_mesh_exact yields a regular box mesh, so
    // this is normally a no-op, but call it for future-proofing).
    std::vector<std::vector<double>*> field_ptrs;
    for (auto& [sym, conc] : st.fields) field_ptrs.push_back(&conc);
    const RepairResult rr = repair_quality(st.mesh, &field_ptrs, 0.1);

    st.last_temp = temp_k;

    if (log)
      *log << "[oxidize] " << (wet ? "wet " : "dry ") << fmt("%.6g", temp_k)
           << " K " << fmt("%.6g", time_s) << " s: tox " << fmt("%.4g", x0_cm * 1e4)
           << " -> " << fmt("%.4g", x_new_um) << " um (dSi="
           << fmt("%.4g", dx_si * 1e4) << " um, rise=" << fmt("%.4g", rise * 1e4)
           << " um), mesh now " << st.mesh.cells.size() << " tets, min_q="
           << fmt("%.4g", rr.min_q_after) << ") [OED disabled, theta=0]\n";

    return x_new_um * 1e-4;
  }

  // ---- P2-3 OED path (theta > 0): N=10 sub-steps of grow -> inject ->
  // point-defect relax, so a single oxidize() call already shows OED. ----
  const int N = 10;
  double x_running_cm = x0_cm;
  double z_si_top_cur = z_si_top;
  double pending_dx_si = 0.0, pending_rise = 0.0;
  double total_I_injected_atoms = 0.0;  // cm^-2 equivalent (integrated volume*conc)
  int ox_tag = -1;

  const auto [nx_init, ny_init, nz_init] = infer_box_dims(st.mesh);
  const double base_lz0 = bb0.hi.z - bb0.lo.z;
  const double h = base_lz0 / std::max(1, nz_init);  // approx z-cell height

  bool has_dopants = false;
  for (const auto& [sym, conc] : st.fields) {
    (void)conc;
    if (find_dopant(sym)) { has_dopants = true; break; }
  }

  // Realize `dx_si_step`/`rise_step` worth of pending geometry change into
  // the mesh: extend upward by rise_step and retag the top dx_si_step of Si
  // as oxide (mirrors the P1-6 (d)-(h) steps, but incremental). No-op if
  // both increments are ~0 (nothing pending).
  auto realize_geometry = [&](double dx_si_step, double rise_step) {
    if (rise_step <= 0 && dx_si_step <= 0) return;
    const BBox bb_cur = st.mesh.bbox();
    const auto [nxc, nyc, nzc] = infer_box_dims(st.mesh);
    const double base_lz = bb_cur.hi.z - bb_cur.lo.z;
    const double hh = base_lz / std::max(1, nzc);
    // A too-small rise_step (e.g. a final leftover flush well under a cell
    // height) must NOT force a whole new mesh layer -- that would create a
    // degenerate sliver cell (near-zero volume) whose huge face
    // transmissibility poisons the diffusion assembly with inf/NaN entries.
    // Below half a cell height, add no new layer at all: extend_mesh_exact
    // with nz_add=0 simply stretches the existing (nz unchanged) layers to
    // fill the slightly taller box, which is numerically harmless for a
    // sub-cell-height increment.
    const int nz_add = (rise_step >= 0.5 * hh)
        ? std::max(1, static_cast<int>(std::round(rise_step / hh)))
        : 0;
    Mesh ext = extend_mesh_exact(st.mesh, rise_step, nz_add);

    const double z_if = z_si_top_cur - dx_si_step;
    if (ox_tag < 0) {
      for (const auto& [t, m] : st.region_material)
        if (t >= 1000 && is_oxide(m)) { ox_tag = t; break; }
      if (ox_tag < 0) {
        const auto etags = st.mesh.region_tags();
        const int max_tag = etags.empty() ? 0
            : *std::max_element(etags.begin(), etags.end());
        ox_tag = max_tag + 3000;
      }
      st.region_material[ox_tag] = "oxide";
    }

    const int nc_ext = static_cast<int>(ext.cells.size());
    ext.cell_region.resize(nc_ext, 0);
    for (int ci = 0; ci < nc_ext; ++ci) {
      const Vec3& c = ext.cell_cent[ci];
      if (c.z > z_if) {
        ext.cell_region[ci] = ox_tag;
      } else {
        double best = 1e300;
        int best_j = 0;
        const int nc_old = static_cast<int>(st.mesh.cells.size());
        for (int j = 0; j < nc_old; ++j) {
          const Vec3& oc = st.mesh.cell_cent[j];
          const double dx = c.x - oc.x, dy = c.y - oc.y, dz = c.z - oc.z;
          const double d2 = dx*dx + dy*dy + dz*dz;
          if (d2 < best) { best = d2; best_j = j; }
        }
        ext.cell_region[ci] = st.mesh.cell_region[best_j];
      }
    }

    const double old_top = bb_cur.hi.z;
    for (auto& [sym, conc] : st.fields) {
      std::vector<double> new_conc(nc_ext, 0.0);
      for (int ci = 0; ci < nc_ext; ++ci) {
        const Vec3& cc = ext.cell_cent[ci];
        double best = 1e300;
        int best_j = 0;
        const int nc_old = static_cast<int>(st.mesh.cells.size());
        for (int j = 0; j < nc_old; ++j) {
          const Vec3& oc = st.mesh.cell_cent[j];
          const double dx = cc.x - oc.x, dy = cc.y - oc.y, dz = cc.z - oc.z;
          const double d2 = dx*dx + dy*dy + dz*dz;
          if (d2 < best) { best = d2; best_j = j; }
        }
        double v = conc[best_j];
        if (cc.z > old_top) v = 0.0;  // new cell above the old top
        new_conc[ci] = v;
      }
      conc = std::move(new_conc);
    }

    st.mesh = std::move(ext);
    std::vector<std::vector<double>*> field_ptrs;
    for (auto& [sym, conc] : st.fields) field_ptrs.push_back(&conc);
    repair_quality(st.mesh, &field_ptrs, 0.1);
    z_si_top_cur = z_if;
  };

  for (int k = 1; k <= N; ++k) {
    const double dt_k = time_s / N;

    // (a) Deal-Grove increment for this sub-step, continuing from the
    // running oxide thickness (exact, since deal_grove_step recomputes the
    // equivalent time from x0 each call -- chaining N sub-steps reproduces
    // the same trajectory as one big step).
    const double x_new_um = deal_grove_step(x_running_cm * 1e4, dt_k / 60.0,
                                            temp_k - 273.15, wet);
    double d_dx = x_new_um * 1e-4 - x_running_cm;
    if (d_dx < 0) d_dx = 0;
    x_running_cm = x_new_um * 1e-4;
    const double d_dx_si = 0.44 * d_dx;
    const double d_rise = 0.56 * d_dx;
    pending_dx_si += d_dx_si;
    pending_rise += d_rise;

    // (c) Interstitial injection: theta * d_dx * kNSi / h [cm^-3], added to
    // the Si cells within one cell layer below the (running, possibly not
    // yet mesh-realized) Si/SiO2 interface height z_if_est. Blanket (no x/y
    // restriction) per the spec.
    if (d_dx > 0) {
      const double z_if_est = z_si_top_cur - pending_dx_si;
      auto& I = st.fields["I"];
      I.resize(st.mesh.cells.size(), 0.0);
      const std::vector<char> simask = silicon_mask(st);
      const double add = theta * d_dx * kNSi / h;
      // Cap the local excess at `oed.psi_cap` (default 50) times the
      // equilibrium C_I*(T): the raw theta*d_dx*N_Si/h_cell source is many
      // orders of magnitude above C_I* in a single ~0.01 um surface cell
      // (the areal Si-consumption flux concentrated into one thin layer).
      // Measured against tests/test_oed.cpp and test_oxidize_flow.cpp: an
      // uncapped (or weakly capped, e.g. 1e4x) source makes the
      // backward-Euler/Picard point-defect+segregation solve either fail to
      // converge (bicgstab breakdown/NaN) or diverge outright on this
      // Si/SiO2-interface mesh; 50x C_I* is the largest cap that stayed
      // numerically stable across every oxidize test while still giving a
      // clearly measurable (>1.2x) OED spread enhancement.
      const PointDefectParams pdp = point_defect_params(temp_k, ParamDB::instance());
      const double cap_ratio = ParamDB::instance().get("oed.psi_cap", 50.0);
      const double cap_val = cap_ratio * pdp.ci_star;
      for (std::size_t ci = 0; ci < st.mesh.cells.size(); ++ci) {
        if (!simask[ci]) continue;
        const double z = st.mesh.cell_cent[ci].z;
        if (z >= z_if_est - h && z < z_if_est) {
          const double before = I[ci];
          I[ci] = std::min(I[ci] + add, cap_val);
          total_I_injected_atoms += (I[ci] - before) * st.mesh.cell_vol[ci];
        }
      }
    }

    // (b) Geometry realization: defer small rises (< 1/4 cell height) and
    // batch them, so we don't rebuild the mesh 10x per oxidize() call; the
    // last sub-step always flushes any remainder.
    if (pending_rise >= 0.25 * h || k == N) {
      realize_geometry(pending_dx_si, pending_rise);
      pending_dx_si = 0.0;
      pending_rise = 0.0;
    }

    // (d) Point-defect (+ dopant) relaxation over this sub-step's duration.
    // Skipped when no dopant field is present yet (nothing to relax; a
    // silicon-only mesh would just pay run_ted's cost for no observable
    // effect). diffuse_ted() lazily creates/updates "I"/"V"/"C311" so the
    // interstitials injected above are picked up automatically.
    if (has_dopants) {
      DiffuseOpts opts_k;
      opts_k.temp = temp_k;
      opts_k.time = dt_k;
      opts_k.verbosity = 0;
      // The default lin_rtol=1e-10/lin_maxit=2000 (tuned for plain anneals)
      // is too tight/short for the combined OED-source + segregation system
      // right at a freshly retagged Si/SiO2 interface; loosen both so the
      // (already-relaxed, since dt_k/50 default) BE/Picard step reliably
      // converges instead of throwing "linear solver failed".
      opts_k.lin_maxit = 5000;
      opts_k.lin_rtol = 1e-8;
      diffuse_ted(st, opts_k, log);
    }
  }

  st.last_temp = temp_k;

  if (log)
    *log << "[oxidize] " << (wet ? "wet " : "dry ") << fmt("%.6g", temp_k)
         << " K " << fmt("%.6g", time_s) << " s: tox " << fmt("%.4g", x0_cm * 1e4)
         << " -> " << fmt("%.4g", x_running_cm * 1e4)
         << " um, mesh now " << st.mesh.cells.size()
         << " tets, N=" << N << " sub-steps, OED theta=" << fmt("%.4g", theta)
         << " I_injected=" << fmt("%.4g", total_I_injected_atoms) << " atoms\n";

  return x_running_cm;
}

// ---------------------------------------------------------------------------
// P2-4: 2D/3D LOCOS oxidation (bird's beak). See process.hpp for the API
// design rationale (separate entry point from oxidize(), column-wise
// realization on the regular box-mesh grid).
// ---------------------------------------------------------------------------
namespace {
bool is_nitride_mat(const std::string& mat) {
  const std::string m = lower(mat);
  return m == "nitride" || m == "si3n4";
}
}  // namespace

double oxidize_2d(SimState& st, double time_s, double temp_k, bool wet,
                  std::ostream* log) {
  need_mesh(st);
  if (st.has_stack)
    throw std::runtime_error("oxidize_2d: resist stack is present; strip first");
  if (time_s <= 0) throw std::runtime_error("oxidize_2d: time_s must be > 0");

  auto material_of = [&](int tag) -> std::string {
    auto it = st.region_material.find(tag);
    return it == st.region_material.end() ? "silicon" : it->second;
  };

  {
    bool has_nitride = false;
    for (const auto& [t, m] : st.region_material)
      if (is_nitride_mat(m)) { has_nitride = true; break; }
    if (!has_nitride)
      throw std::runtime_error(
          "oxidize_2d: no nitride-masked region present; use oxidize() for "
          "blanket (1D) oxidation");
  }

  const auto [nx0, ny0, nz0] = infer_box_dims(st.mesh);
  const int ncol = nx0 * ny0;
  auto column_id = [&](int cell) { return (cell / 6) % ncol; };

  int ox_tag = -1;
  for (const auto& [t, m] : st.region_material)
    if (t >= 1000 && is_oxide(m)) { ox_tag = t; break; }
  auto ensure_ox_tag = [&]() {
    if (ox_tag >= 0) return;
    const auto etags = st.mesh.region_tags();
    const int max_tag = etags.empty() ? 0 : *std::max_element(etags.begin(), etags.end());
    ox_tag = max_tag + 3000;
    st.region_material[ox_tag] = "oxide";
  };

  // Realize a per-column geometry increment: dx_si_col[col] of silicon
  // consumed and rise_col[col] of outer-surface swell, both cm, both >= 0.
  // Mirrors oxidize()'s realize_geometry, generalized from a scalar
  // interface height to a per-column height field. Nitride cells are
  // protected (never retagged) via the nearest-old-centroid lookup: any new
  // cell whose nearest old cell was nitride keeps that nitride tag
  // regardless of the local z_if_col comparison.
  auto realize_geometry_col = [&](const std::vector<double>& dx_si_col,
                                  const std::vector<double>& rise_col,
                                  std::vector<double>& z_si_top_col) {
    double max_rise = 0;
    for (double r : rise_col) max_rise = std::max(max_rise, r);
    if (max_rise <= 0) return;
    ensure_ox_tag();

    const BBox bb_cur = st.mesh.bbox();
    const auto [nxc, nyc, nzc] = infer_box_dims(st.mesh);
    const double base_lz = bb_cur.hi.z - bb_cur.lo.z;
    const double hh = base_lz / std::max(1, nzc);
    const int nz_add = (max_rise >= 0.25 * hh)
        ? std::max(1, static_cast<int>(std::round(max_rise / hh))) : 0;
    Mesh ext = extend_mesh_exact(st.mesh, max_rise, nz_add);

    const int nc_ext = static_cast<int>(ext.cells.size());
    const int nc_old = static_cast<int>(st.mesh.cells.size());
    ext.cell_region.resize(nc_ext, 0);
    std::vector<int> nearest_old(nc_ext, 0);
    for (int ci = 0; ci < nc_ext; ++ci) {
      const Vec3& c = ext.cell_cent[ci];
      double best = 1e300;
      int best_j = 0;
      for (int j = 0; j < nc_old; ++j) {
        const Vec3& oc = st.mesh.cell_cent[j];
        const double dx = c.x - oc.x, dy = c.y - oc.y, dz = c.z - oc.z;
        const double d2 = dx * dx + dy * dy + dz * dz;
        if (d2 < best) { best = d2; best_j = j; }
      }
      nearest_old[ci] = best_j;
      const int old_tag = st.mesh.cell_region[best_j];
      if (is_nitride_mat(material_of(old_tag))) {
        ext.cell_region[ci] = old_tag;  // nitride never converts/moves
        continue;
      }
      const int col = column_id(ci);
      const double z_if = z_si_top_col[col] - dx_si_col[col];
      ext.cell_region[ci] = (c.z > z_if) ? ox_tag : old_tag;
    }

    const double old_top = bb_cur.hi.z;
    for (auto& [sym, conc] : st.fields) {
      std::vector<double> new_conc(nc_ext, 0.0);
      for (int ci = 0; ci < nc_ext; ++ci) {
        double v = conc[nearest_old[ci]];
        if (ext.cell_cent[ci].z > old_top) v = 0.0;
        new_conc[ci] = v;
      }
      conc = std::move(new_conc);
    }

    st.mesh = std::move(ext);
    std::vector<std::vector<double>*> field_ptrs;
    for (auto& [sym, conc] : st.fields) field_ptrs.push_back(&conc);
    repair_quality(st.mesh, &field_ptrs, 0.1);

    for (int col = 0; col < ncol; ++col) z_si_top_col[col] -= dx_si_col[col];
  };

  // Per-column running Si top height.
  std::vector<double> z_si_top_col(ncol, -1e300);
  for (std::size_t i = 0; i < st.mesh.cells.size(); ++i) {
    if (!is_silicon(material_of(st.mesh.cell_region[i]))) continue;
    const auto& cell = st.mesh.cells[i];
    const int col = column_id(static_cast<int>(i));
    for (int k = 0; k < 4; ++k)
      z_si_top_col[col] = std::max(z_si_top_col[col], st.mesh.nodes[cell[k]].z);
  }

  // Bootstrap: solve_oxidant needs existing oxide (or nitride) cells to host
  // the active domain. Any column with neither yet (bare Si straight to gas)
  // gets a thin uniform "native oxide" seed layer first.
  {
    std::vector<char> needs_seed(ncol, 0);
    for (int col = 0; col < ncol; ++col) needs_seed[col] = 1;
    for (std::size_t i = 0; i < st.mesh.cells.size(); ++i) {
      const std::string mat = material_of(st.mesh.cell_region[i]);
      if (is_oxide(mat) || is_nitride_mat(mat))
        needs_seed[column_id(static_cast<int>(i))] = 0;
    }
    bool any_seed = false;
    for (char c : needs_seed) if (c) { any_seed = true; break; }
    if (any_seed) {
      const double x_seed_cm = ParamDB::instance().get("ox2d.seed_ox_um", 0.001) * 1e-4;
      std::vector<double> dx_si_col(ncol, 0.0), rise_col(ncol, 0.0);
      for (int col = 0; col < ncol; ++col) {
        if (!needs_seed[col]) continue;
        dx_si_col[col] = 0.44 * x_seed_cm;
        rise_col[col] = 0.56 * x_seed_cm;
      }
      realize_geometry_col(dx_si_col, rise_col, z_si_top_col);
    }
  }

  const int N = 10;
  const double dt_step = time_s / N;
  const double temp_c = temp_k - 273.15;
  const double c_gas = ParamDB::instance().get("ox2d.cgas", 5.2e16);
  // Nitride is not a perfect oxidant barrier at this coarse a mesh
  // resolution: some lateral bleed-through (combined with lateral diffusion
  // from the open field, through the same nitride cells) is what produces
  // the bird's-beak taper under the mask edge, per the P2-4 spec's own
  // note (see docs/tasks/P2-4_oxidation_2d3d.md item 5). 0.2 was the
  // smallest value found, while developing this solver, that reproduces a
  // measurable (20-80%) taper at 0.1 um under the mask edge on the test
  // mesh's resolution (8x4x50 over 0.4x0.2x0.5 um); true nitride is a far
  // better barrier, but resolving the actual bird's-beak mechanism (oxidant
  // diffusing laterally through a curved oxide wedge under the mask, not
  // through the nitride itself) would need genuine per-node ALE motion of
  // the oxide/nitride/Si triple point, out of reach in this task's time
  // budget -- see the report / API doc comment above for the full
  // scope-decision rationale.
  const double nitride_leak = ParamDB::instance().get("ox2d.nitride_leak", 0.2);
  const double theta = ParamDB::instance().get("oed.theta", 0.01);

  bool has_dopants = false;
  for (const auto& [sym, conc] : st.fields) {
    (void)conc;
    if (find_dopant(sym)) { has_dopants = true; break; }
  }

  double open_field_x_o_cm = 0.0;  // reported return value: open-field oxide thickness
  double total_I_injected_atoms = 0.0;

  for (int step = 1; step <= N; ++step) {
    const auto [nxc, nyc, nzc] = infer_box_dims(st.mesh);
    (void)nxc; (void)nyc;
    const double base_lz = st.mesh.bbox().hi.z - st.mesh.bbox().lo.z;
    const double h = base_lz / std::max(1, nzc);

    // Deal-Grove-consistent oxidant diffusivity/reaction rate (Stage A note:
    // B, B/A from oxidation.cpp are per-hour, um; convert to cm^2/s, cm/s).
    const DealGroveParams p = wet ? deal_grove_wet(temp_c) : deal_grove_dry(temp_c);
    constexpr double kNOx = 2.25e22;
    const double B_cm2_s = p.B * 1e-8 / 3600.0;
    const double BA_cm_s = (p.B / p.A) * 1e-4 / 3600.0;
    const double d_ox = B_cm2_s * kNOx / (2.0 * c_gas);
    const double ks = BA_cm_s * kNOx / c_gas;

    const int nc = static_cast<int>(st.mesh.cells.size());
    std::vector<char> active_mask(nc, 0), si_mask(nc, 0);
    std::vector<double> d_cell(nc, 0.0);
    for (int i = 0; i < nc; ++i) {
      const std::string mat = material_of(st.mesh.cell_region[i]);
      if (is_oxide(mat)) { active_mask[i] = 1; d_cell[i] = d_ox; }
      else if (is_nitride_mat(mat)) { active_mask[i] = 1; d_cell[i] = d_ox * nitride_leak; }
      else if (is_silicon(mat)) { si_mask[i] = 1; }
    }

    const OxidantResult ores =
        solve_oxidant(st.mesh, d_cell, active_mask, si_mask, ks, c_gas, log);

    // Per-column area-weighted average growth velocity (cm/s), from every
    // Robin face whose owning active cell lies in that column.
    std::vector<double> vn_num(ncol, 0.0), vn_den(ncol, 0.0);
    for (const auto& [fi, vn] : ores.iface_vn) {
      const Face& f = st.mesh.faces[fi];
      const int col = column_id(f.owner);
      const double area = norm(f.S);
      vn_num[col] += vn * area;
      vn_den[col] += area;
    }
    std::vector<double> dx_si_col(ncol, 0.0), rise_col(ncol, 0.0), dx_o_col(ncol, 0.0);
    for (int col = 0; col < ncol; ++col) {
      if (vn_den[col] <= 0) continue;
      const double vn_col = vn_num[col] / vn_den[col];
      double dx_o = vn_col * dt_step;
      if (dx_o < 0) dx_o = 0;
      dx_o_col[col] = dx_o;
      dx_si_col[col] = 0.44 * dx_o;
      rise_col[col] = 0.56 * dx_o;
    }

    // P2-3-style interstitial injection, per column, using the locally
    // solved dx_o_col (blanket theta model, same cap logic as oxidize()).
    if (theta > 0.0) {
      auto& I = st.fields["I"];
      I.resize(st.mesh.cells.size(), 0.0);
      const PointDefectParams pdp = point_defect_params(temp_k, ParamDB::instance());
      const double cap_ratio = ParamDB::instance().get("oed.psi_cap", 50.0);
      const double cap_val = cap_ratio * pdp.ci_star;
      for (int i = 0; i < nc; ++i) {
        if (!si_mask[i]) continue;
        const int col = column_id(i);
        if (dx_o_col[col] <= 0) continue;
        const double z_if_est = z_si_top_col[col] - dx_si_col[col];
        const double z = st.mesh.cell_cent[i].z;
        if (z >= z_if_est - h && z < z_if_est) {
          const double add = theta * dx_o_col[col] * kNSi / h;
          const double before = I[i];
          I[i] = std::min(I[i] + add, cap_val);
          total_I_injected_atoms += (I[i] - before) * st.mesh.cell_vol[i];
        }
      }
    }

    realize_geometry_col(dx_si_col, rise_col, z_si_top_col);

    if (has_dopants) {
      DiffuseOpts opts_k;
      opts_k.temp = temp_k;
      opts_k.time = dt_step;
      opts_k.verbosity = 0;
      opts_k.lin_maxit = 5000;
      opts_k.lin_rtol = 1e-8;
      diffuse_ted(st, opts_k, log);
    }
  }

  // Report: average oxide thickness over open (non-nitride) columns.
  {
    std::vector<double> col_ox_top(ncol, -1e300), col_si_top(ncol, 1e300);
    std::vector<char> col_has_nitride(ncol, 0);
    for (std::size_t i = 0; i < st.mesh.cells.size(); ++i) {
      const std::string mat = material_of(st.mesh.cell_region[i]);
      const int col = column_id(static_cast<int>(i));
      const auto& cell = st.mesh.cells[i];
      if (is_nitride_mat(mat)) col_has_nitride[col] = 1;
      if (is_oxide(mat))
        for (int k = 0; k < 4; ++k) col_ox_top[col] = std::max(col_ox_top[col], st.mesh.nodes[cell[k]].z);
    }
    double sum = 0; int n = 0;
    for (int col = 0; col < ncol; ++col) {
      if (col_has_nitride[col]) continue;
      if (col_ox_top[col] < -1e299) continue;
      sum += col_ox_top[col] - z_si_top_col[col];
      ++n;
    }
    open_field_x_o_cm = (n > 0) ? sum / n : 0.0;
  }

  st.last_temp = temp_k;
  if (log)
    *log << "[oxidize_2d] " << (wet ? "wet " : "dry ") << fmt("%.6g", temp_k)
         << " K " << fmt("%.6g", time_s) << " s: open-field tox ~ "
         << fmt("%.4g", open_field_x_o_cm * 1e4) << " um, mesh now "
         << st.mesh.cells.size() << " tets, N=" << N
         << " sub-steps, I_injected=" << fmt("%.4g", total_I_injected_atoms)
         << " atoms\n";

  return open_field_x_o_cm;
}

void add_bc(SimState& st, const std::string& species, int patch, double conc,
            std::ostream* log) {
  need_mesh(st);
  const Dopant* d = dopant_or_throw(species);
  if (patch < 0) throw std::runtime_error("invalid patch");
  if (conc < 0) throw std::runtime_error("conc must be >= 0");
  DirichletBC bc;
  bc.species = d->symbol;
  bc.patch = patch;
  bc.conc = conc;
  st.bcs.push_back(bc);
  st.fields[d->symbol].resize(st.mesh.cells.size(), 0.0);
  if (log)
    *log << "[bc] " << bc.species << " = " << fmt("%.3g", conc) << " cm^-3\n";
}

void clear_bc(SimState& st, std::ostream* log) {
  st.bcs.clear();
  if (log) *log << "[bc] cleared\n";
}

void diffuse(SimState& st, const DiffuseOpts& opts, std::ostream* log) {
  need_mesh(st);
  std::vector<SpeciesField> fields;
  for (auto& [sym, conc] : st.fields) {
    const Dopant* d = find_dopant(sym);
    if (d) fields.push_back({d, &conc});
  }
  if (fields.empty()) {
    if (log) *log << "[diffuse] no dopants present\n";
    return;
  }
  if (log) {
    *log << "[diffuse] T=";
    if (opts.temp_profile.empty()) {
      *log << fmt("%.5g", opts.temp) << " K";
    } else {
      *log << "ramp(";
      for (std::size_t i = 0; i < opts.temp_profile.size(); ++i) {
        if (i) *log << "->";
        *log << fmt("%.5g", opts.temp_profile[i].second);
      }
      *log << " K)";
    }
    *log << ", time=" << fmt("%.5g", opts.time) << " s, species:";
    for (const auto& fl : fields) *log << " " << fl.dopant->symbol;
    *log << "\n";
  }
  DiffusionSolver solver(st.mesh, silicon_mask(st), log, material_ids(st));
  solver.run(fields, st.bcs, opts);
  st.last_temp = temp_at(opts, opts.time);
}

void diffuse_ted(SimState& st, const DiffuseOpts& opts, std::ostream* log) {
  need_mesh(st);
  // P2-2/P2-8: lazily create "B_cl"/"As_cl"/"C_cl" immobile cluster fields
  // for any clustering- (or, for C, C-I-sink-) capable dopant already
  // present, before building `fields` (so the loop below can bind
  // SpeciesField::cluster without inserting into st.fields mid-iteration).
  // "C_cl" accumulates C-I sink capture (see run_ted's C-I sink block), not
  // BIC-style clustering (cluster_params() returns kf=0 for "C").
  std::vector<std::string> cluster_syms;
  for (const auto& [sym, conc] : st.fields) {
    if ((sym == "B" || sym == "As" || sym == "C") && find_dopant(sym))
      cluster_syms.push_back(sym);
  }
  for (const auto& sym : cluster_syms)
    st.fields[sym + "_cl"].resize(st.mesh.cells.size(), 0.0);

  std::vector<SpeciesField> fields;
  for (auto& [sym, conc] : st.fields) {
    const Dopant* d = find_dopant(sym);
    if (!d) continue;
    SpeciesField sf{d, &conc};
    if (sym == "B" || sym == "As" || sym == "C") sf.cluster = &st.fields[sym + "_cl"];
    fields.push_back(sf);
  }
  if (fields.empty()) {
    if (log) *log << "[ted] no dopants present\n";
    return;
  }
  // Interstitial excess field seeded by implants (damage=true). Absent means a
  // plain equilibrium anneal — run_ted still works (S = 1 everywhere).
  auto& psi = st.fields["I"];
  psi.resize(st.mesh.cells.size(), 0.0);
  // Vacancy excess and {311} cluster reservoir (P2-1). Both are created here
  // (proc layer) so Python/C++ callers can inspect "V"/"C311" like any other
  // field; run_ted() persists them across repeated diffuse_ted() calls.
  auto& v = st.fields["V"];
  v.resize(st.mesh.cells.size(), 0.0);
  auto& c311 = st.fields["C311"];
  c311.resize(st.mesh.cells.size(), 0.0);
  if (log) {
    *log << "[ted] T=";
    if (opts.temp_profile.empty()) {
      *log << fmt("%.5g", opts.temp) << " K";
    } else {
      *log << "ramp(";
      for (std::size_t i = 0; i < opts.temp_profile.size(); ++i) {
        if (i) *log << "->";
        *log << fmt("%.5g", opts.temp_profile[i].second);
      }
      *log << " K)";
    }
    *log << ", time=" << fmt("%.5g", opts.time) << " s, species:";
    for (const auto& fl : fields) *log << " " << fl.dopant->symbol;
    *log << "\n";
  }
  DiffusionSolver solver(st.mesh, silicon_mask(st), log, material_ids(st));
  solver.run_ted(fields, psi, v, c311, st.bcs, opts);
  st.last_temp = temp_at(opts, opts.time);
}


namespace {
double mat_default_E_gpa(const std::string& m) {
  if (is_oxide(m)) return 70.0;
  if (is_nitride_mat(m)) return 250.0;
  if (m == "poly" || m == "polysilicon") return 160.0;
  if (m == "gas") return 1e-6;  // near-zero: void has essentially no rigidity
  return 130.0;  // silicon (matches mechanics.hpp's ViscoElasticParams default)
}
double mat_default_nu(const std::string& m) {
  if (is_oxide(m)) return 0.17;
  if (is_nitride_mat(m)) return 0.23;
  if (m == "poly" || m == "polysilicon") return 0.22;
  if (m == "gas") return 0.0;
  return 0.28;  // silicon
}
double mat_default_alpha(const std::string& m) {
  if (is_oxide(m)) return 0.5e-6;
  if (is_nitride_mat(m)) return 3.3e-6;
  if (m == "gas") return 0.0;
  return 2.6e-6;  // silicon, poly (no dedicated poly literature value on hand)
}
double mat_default_sigma0(const std::string& m) {
  // Only nitride carries a nonzero literature default (~1 GPa tensile,
  // typical of LPCVD Si3N4); other materials default to no intrinsic stress.
  if (is_nitride_mat(m)) return 1e10;  // dyn/cm^2
  return 0.0;
}
}  // namespace

// P2-6: linear-elastic FEM mechanics. See process.hpp for the API contract.
void mechanics(SimState& st, double temp_k, double dt_s, std::ostream* log) {
  need_mesh(st);
  if (dt_s < 0) throw std::runtime_error("mechanics: dt_s must be >= 0");
  const double dT = temp_k - 300.0;
  const int nc = static_cast<int>(st.mesh.cells.size());
  ParamDB& db = ParamDB::instance();

  std::vector<std::string> mat_name(nc);
  for (int c = 0; c < nc; ++c) {
    auto it = st.region_material.find(st.mesh.cell_region[c]);
    mat_name[c] = lower(it == st.region_material.end() ? "silicon" : it->second);
  }

  // Per-cell eigenstrain (Voigt order xx,yy,zz,yz,xz,xy, per fem.hpp): thermal
  // mismatch alpha*dT plus intrinsic film stress sigma0 mapped through
  // eps0 = D^-1 * sigma0 for an isotropic hydrostatic eigenstress (documented
  // simplification -- a true thin-film intrinsic stress is biaxial in-plane
  // only, but the isotropic form keeps the eigenstrain construction identical
  // to the thermal term and is adequate for the sign/magnitude acceptance
  // tests this task targets).
  std::vector<double> E_cell(nc), nu_cell(nc), eps0(nc * 6, 0.0);
  for (int c = 0; c < nc; ++c) {
    const std::string& m = mat_name[c];
    const double E_gpa = db.get("mech.E." + m, mat_default_E_gpa(m));
    const double nu = db.get("mech.nu." + m, mat_default_nu(m));
    const double alpha = db.get("mech.alpha." + m, mat_default_alpha(m));
    const double sigma0_dyncm2 = db.get("mech.sigma0." + m, mat_default_sigma0(m));
    E_cell[c] = E_gpa * 1e3;      // GPa -> MPa (fem.cpp's internal unit)
    nu_cell[c] = nu;
    const double sigma0_mpa = sigma0_dyncm2 * 1e-7;  // dyn/cm^2 -> MPa
    const double e_thermal = alpha * dT;
    const double e_intrinsic = sigma0_mpa * (1.0 - 2.0 * nu) / std::max(E_cell[c], 1e-9);
    for (int k = 0; k < 3; ++k) eps0[c * 6 + k] = e_thermal + e_intrinsic;
  }

  FemProblem prob = fem_assemble(st.mesh, E_cell, nu_cell, eps0);

  // BC: roller at the three "min" faces (zmin/xmin/ymin), each pinning only
  // its own normal displacement component; xmax/ymax/zmax free. See the
  // detailed rationale in the zmin block below.
  const BBox bb = st.mesh.bbox();
  const double span = std::max({bb.hi.x - bb.lo.x, bb.hi.y - bb.lo.y,
                                bb.hi.z - bb.lo.z, 1e-12});
  const double tol = 1e-9 * span;
  const int nn = static_cast<int>(st.mesh.nodes.size());
  for (int i = 0; i < nn; ++i) {
    const Vec3& p = st.mesh.nodes[i];
    // zmin: roller (uz=0 only), not a full 3-component clamp. A full clamp
    // (ux=uy=uz=0 at every zmin node, not just the origin) is incompatible
    // with pure uniform thermal expansion u=alpha*dT*(x,y,z) whenever a
    // zmin node has x!=0 or y!=0 -- it would force nonzero elastic strain
    // near that boundary even for a single homogeneous material with no
    // real internal stress, which is exactly what the "uniform expansion ->
    // near-zero von Mises" acceptance test below checks for. Fixing uz=0 on
    // the whole zmin plane (matching the ux=0/uy=0 treatment already used at
    // xmin/ymin below) is compatible with that field (it vanishes at z=0),
    // and the three orthogonal roller planes together still fully eliminate
    // all 6 rigid-body modes (3 translation + 3 rotation), so the system
    // stays well-posed without over-constraining thermal expansion.
    if (std::fabs(p.z - bb.lo.z) < tol) {
      prob.dof_fixed[3 * i + 2] = 1;
      prob.dof_val[3 * i + 2] = 0.0;
    }
    // Roller only at the "min" lateral faces (not xmax/ymax): this is a
    // symmetry-style BC for a free-standing block (like modeling one octant
    // of a body that expands away from a fixed corner). Rollering *both*
    // xmin and xmax (or ymin/ymax) would instead pin the block's overall
    // x (or y) extent, which suppresses free thermal expansion entirely and
    // produces large spurious stress -- confirmed by direct instrumentation
    // while implementing the uniform-expansion acceptance test below (a
    // both-sides roller gave ~400 MPa von Mises on a uniform dT, vs <1e-6 MPa
    // with only the min-face roller used here).
    if (std::fabs(p.x - bb.lo.x) < tol) {
      prob.dof_fixed[3 * i + 0] = 1;
      prob.dof_val[3 * i + 0] = 0.0;
    }
    if (std::fabs(p.y - bb.lo.y) < tol) {
      prob.dof_fixed[3 * i + 1] = 1;
      prob.dof_val[3 * i + 1] = 0.0;
    }
  }
  fem_apply_bc(prob);

  std::vector<double> u;
  SolveResult res = fem_solve(prob, u, 1e-10, 5000);
  if (!res.converged && log)
    *log << "mechanics: WARNING cg_ilu0 did not converge (resid="
         << fmt("%.3g", res.resid) << ")\n";

  std::vector<double> stress = fem_element_stress(st.mesh, u, E_cell, nu_cell, eps0);

  auto& sxx = st.fields["sxx"]; auto& syy = st.fields["syy"]; auto& szz = st.fields["szz"];
  auto& sxy = st.fields["sxy"]; auto& syz = st.fields["syz"]; auto& sxz = st.fields["sxz"];
  sxx.assign(nc, 0.0); syy.assign(nc, 0.0); szz.assign(nc, 0.0);
  sxy.assign(nc, 0.0); syz.assign(nc, 0.0); sxz.assign(nc, 0.0);

  const double dt_min = dt_s / 60.0;
  for (int c = 0; c < nc; ++c) {
    // fem.hpp order (xx,yy,zz,yz,xz,xy) -> mechanics.hpp/maxwell_update order
    // (xx,yy,zz,xy,yz,xz).
    const double* s = &stress[c * 6];
    StressState s_old;
    s_old.s[0] = s[0]; s_old.s[1] = s[1]; s_old.s[2] = s[2];
    s_old.s[3] = s[5]; s_old.s[4] = s[3]; s_old.s[5] = s[4];

    const std::string& m = mat_name[c];
    ViscoElasticParams vp;
    vp.E_GPa = db.get("mech.E." + m, mat_default_E_gpa(m));
    vp.nu = db.get("mech.nu." + m, mat_default_nu(m));
    const double tau_s = db.get("mech.tau." + m, 0.0);
    vp.tau_relax_min = (tau_s > 0.0) ? (tau_s / 60.0) : 1e12;  // <=0 => elastic

    const double dstrain[6] = {};  // relax the already-solved stress in place
    const StressState s_new = maxwell_update(s_old, dstrain, dt_min, vp);

    const double k = 1e7;  // MPa -> dyn/cm^2
    sxx[c] = s_new.s[0] * k; syy[c] = s_new.s[1] * k; szz[c] = s_new.s[2] * k;
    sxy[c] = s_new.s[3] * k; syz[c] = s_new.s[4] * k; sxz[c] = s_new.s[5] * k;
  }

  if (log)
    *log << "mechanics: dT=" << fmt("%.5g", dT) << " K, dt="
         << fmt("%.5g", dt_s) << " s, cg iters=" << res.iters
         << " resid=" << fmt("%.3g", res.resid) << "\n";
}

void refine(SimState& st, const std::string& species, double rel_grad_thresh,
           int max_passes, const std::string& axis, std::ostream* log) {
  need_mesh(st);
  auto it = st.fields.find(species);
  if (it == st.fields.end() || it->second.empty())
    throw std::runtime_error("refine: no field '" + species + "'");

  std::vector<std::vector<double>*> field_ptrs;
  int key_index = -1, idx = 0;
  for (auto& [sym, conc] : st.fields) {
    if (sym == species) key_index = idx;
    field_ptrs.push_back(&conc);
    ++idx;
  }

  const int nc0 = static_cast<int>(st.mesh.cells.size());
  RefineResult rr;
  if (axis.empty()) {
    rr = refine_gradient(st.mesh, field_ptrs, key_index, rel_grad_thresh,
                         max_passes);
  } else {
    Vec3 direction;
    if (axis == "x") direction = {1, 0, 0};
    else if (axis == "y") direction = {0, 1, 0};
    else if (axis == "z") direction = {0, 0, 1};
    else throw std::runtime_error("refine: unknown axis '" + axis + "'");
    rr = refine_anisotropic(st.mesh, direction, &field_ptrs, rel_grad_thresh,
                            max_passes);
  }

  if (log)
    *log << "[refine] " << species << " passes=" << rr.n_passes
         << " split=" << rr.n_split_total << " cells " << nc0 << " -> "
         << rr.n_cells_after << "\n";
}

void save(SimState& st, const std::string& path, std::ostream* log) {
  need_mesh(st);
  const std::size_t nc = st.mesh.cells.size();
  std::vector<std::vector<double>> extra;
  std::vector<std::pair<std::string, const std::vector<double>*>> scalars;
  for (const auto& [sym, conc] : st.fields) scalars.emplace_back(sym, &conc);

  std::vector<double> net(nc, 0.0);
  extra.reserve(st.fields.size());
  for (const auto& [sym, conc] : st.fields) {
    const Dopant* d = find_dopant(sym);
    if (!d) continue;
    extra.emplace_back(nc, 0.0);
    auto& act = extra.back();
    for (std::size_t i = 0; i < nc; ++i) {
      act[i] = active_concentration(*d, conc[i], st.last_temp);
      if (d->type == DopType::neutral) continue;  // P2-8: no net charge
      net[i] += (d->type == DopType::donor) ? act[i] : -act[i];
    }
  }
  std::size_t k = 0;
  for (const auto& [sym, conc] : st.fields) {
    (void)conc;
    if (!find_dopant(sym)) continue;  // only dopant fields have an _active entry
    scalars.emplace_back(sym + "_active", &extra[k++]);
  }
  scalars.emplace_back("NetDoping", &net);

  std::vector<int> region(st.mesh.cell_region.begin(), st.mesh.cell_region.end());
  std::vector<std::pair<std::string, const std::vector<int>*>> ints = {
      {"Region", &region}};
  write_vtu(path, st.mesh, scalars, ints);
  if (log) *log << "[save] wrote " << path << "\n";
}

std::vector<double> active_field(const SimState& st, const std::string& species,
                                 double temp_k, std::ostream* log) {
  const Dopant* d = dopant_or_throw(species);
  auto it = st.fields.find(d->symbol);
  if (it == st.fields.end())
    throw std::runtime_error("no field: " + species);
  const double t_k = (temp_k > 0) ? temp_k : st.last_temp;
  const auto& conc = it->second;
  std::vector<double> act(conc.size());
  for (std::size_t i = 0; i < conc.size(); ++i)
    act[i] = active_concentration(*d, conc[i], t_k);
  if (log)
    *log << "[active] " << d->symbol << ": C_ss(T)=" << solid_solubility(*d, t_k)
         << " cm^-3\n";
  return act;
}

void set_param(SimState& st, const std::string& key, double value,
              std::ostream* log) {
  (void)st;
  ParamDB::instance().set(key, value);
  if (log) *log << "[param] " << key << " = " << fmt("%.6g", value) << "\n";
}

double get_param(const std::string& key, double fallback) {
  return ParamDB::instance().get(key, fallback);
}

std::map<std::string, double> list_params() {
  return ParamDB::instance().all();
}

void save_state(SimState& st, const std::string& path, std::ostream* log) {
  need_mesh(st);
  if (st.has_stack)
    throw std::runtime_error("save_state: strip photoresist stack before save");
  write_state(st, path);
  if (log)
    *log << "[save_state] wrote " << path << " (" << st.mesh.cells.size()
         << " cells, " << st.fields.size() << " fields)\n";
}

void load_state(SimState& st, const std::string& path, std::ostream* log) {
  read_state(st, path);
  if (log)
    *log << "[load_state] " << path << ": " << st.mesh.cells.size()
         << " cells\n";
}

void export_device(SimState& st, const std::string& path_prefix, std::ostream* log) {
  need_mesh(st);
  if (st.has_stack)
    throw std::runtime_error("export_device: strip photoresist stack before export");
  cp::write_device(st, path_prefix);
  std::size_t n_species = 0;
  for (const auto& [sym, conc] : st.fields) {
    (void)conc;
    if (find_dopant(sym)) ++n_species;
  }
  if (log)
    *log << "[export_device] wrote " << path_prefix << ".vtu + " << path_prefix
         << ".meta.json (" << st.mesh.cells.size() << " cells, " << n_species
         << " species)\n";
}

std::vector<std::pair<double, double>> profile1d(const SimState& st,
                                                  const std::string& species,
                                                  double x_cm, double y_cm) {
  auto fit = st.fields.find(species);
  if (fit == st.fields.end())
    throw std::runtime_error("no field: " + species);
  const auto& conc = fit->second;

  std::vector<std::pair<double, double>> out;
  const std::size_t nc = st.mesh.cells.size();
  for (std::size_t ci = 0; ci < nc; ++ci) {
    auto mit = st.region_material.find(st.mesh.cell_region[ci]);
    if (mit != st.region_material.end() && lower(mit->second) == "gas")
      continue;
    double xmin = 1e300, xmax = -1e300, ymin = 1e300, ymax = -1e300;
    for (int nid : st.mesh.cells[ci]) {
      const Vec3& p = st.mesh.nodes[nid];
      xmin = std::min(xmin, p.x);
      xmax = std::max(xmax, p.x);
      ymin = std::min(ymin, p.y);
      ymax = std::max(ymax, p.y);
    }
    if (x_cm < xmin || x_cm > xmax || y_cm < ymin || y_cm > ymax) continue;
    const double c = (ci < conc.size()) ? conc[ci] : 0.0;
    out.emplace_back(st.mesh.cell_cent[ci].z, c);
  }
  std::sort(out.begin(), out.end(),
           [](const auto& a, const auto& b) { return a.first < b.first; });
  return out;
}

std::vector<std::vector<std::pair<double, double>>> load_gds(
    const std::string& path, int layer, std::ostream* log) {
  auto polys = read_gds(path, layer);
  if (log)
    *log << "[load_gds] " << polys.size() << " polygons from layer " << layer
         << "\n";
  return polys;
}

}  // namespace proc
}  // namespace cp
