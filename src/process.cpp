#include "cprocess/process.hpp"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <map>
#include <ostream>
#include <stdexcept>
#include <utility>

#include "cprocess/field_transfer.hpp"
#include "cprocess/gds_reader.hpp"
#include "cprocess/gmsh_reader.hpp"
#include "cprocess/oxidation.hpp"
#include "cprocess/param_db.hpp"
#include "cprocess/remesh.hpp"
#include "cprocess/state_io.hpp"
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

  // (c) Deal-Grove increment (unit conversion: cm -> um, s -> min, K -> C).
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

  // (h) Quality repair (extend_mesh_exact yields a regular box mesh, so this
  // is normally a no-op, but call it for future-proofing).
  std::vector<std::vector<double>*> field_ptrs;
  for (auto& [sym, conc] : st.fields) field_ptrs.push_back(&conc);
  const RepairResult rr = repair_quality(st.mesh, &field_ptrs, 0.1);

  // (i)
  st.last_temp = temp_k;

  if (log)
    *log << "[oxidize] " << (wet ? "wet " : "dry ") << fmt("%.6g", temp_k)
         << " K " << fmt("%.6g", time_s) << " s: tox " << fmt("%.4g", x0_cm * 1e4)
         << " -> " << fmt("%.4g", x_new_um) << " um (dSi="
         << fmt("%.4g", dx_si * 1e4) << " um, rise=" << fmt("%.4g", rise * 1e4)
         << " um), mesh now " << st.mesh.cells.size() << " tets, min_q="
         << fmt("%.4g", rr.min_q_after) << "\n";

  return x_new_um * 1e-4;
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

void refine(SimState& st, const std::string& species, double rel_grad_thresh,
           int max_passes, std::ostream* log) {
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
  RefineResult rr = refine_gradient(st.mesh, field_ptrs, key_index,
                                    rel_grad_thresh, max_passes);

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
