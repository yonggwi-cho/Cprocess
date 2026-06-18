#include "cprocess/process.hpp"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <ostream>
#include <stdexcept>
#include <utility>

#include "cprocess/field_transfer.hpp"
#include "cprocess/gmsh_reader.hpp"
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

const Dopant* dopant_or_throw(const std::string& species) {
  const Dopant* d = find_dopant(species);
  if (!d) throw std::runtime_error("unknown species '" + species + "' (B, P, As, Sb)");
  return d;
}

void need_mesh(const SimState& st) {
  if (!st.has_mesh) throw std::runtime_error("no mesh defined yet");
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
                     std::ostream* log) {
  need_mesh(st);
  ImplantParams p;
  p.dopant = dopant_or_throw(species);
  p.dose = dose;
  if (energy_kev > 0) {
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
  const double atoms = apply_implant(st.mesh, silicon_mask(st), p, f);
  if (log)
    *log << "[implant] " << p.dopant->symbol << " gauss: dose="
         << fmt("%.3g", p.dose) << " cm^-2, Rp=" << fmt("%.4g", p.rp * 1e4)
         << " um, dRp=" << fmt("%.4g", p.drp * 1e4) << " um\n";
  return atoms;
}

McImplantStats implant_mc(SimState& st, const std::string& species, double dose,
                          double energy_kev, long long ions, double tilt_deg,
                          double rotation_deg, unsigned long long seed,
                          int threads, bool channeling, bool has_window,
                          double x1, double x2, double y1, double y2,
                          std::ostream* log) {
  need_mesh(st);
  const Dopant* dop = dopant_or_throw(species);
  auto& f = st.fields[dop->symbol];
  f.resize(st.mesh.cells.size(), 0.0);

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
    const McImplantStats s =
        apply_mc_implant(st.stack, stack_si, p, stack_conc, nullptr);

    const std::vector<double> transferred =
        transfer_field_nearest(st.stack, stack_conc, st.mesh);
    for (std::size_t i = 0; i < f.size(); ++i) f[i] += transferred[i];

    if (log) {
      const double n = static_cast<double>(p.ions);
      *log << "[implant] " << dop->symbol << " MC (physical resist): E="
           << fmt("%.4g", p.energy_kev) << " keV, deposited_in_Si="
           << fmt("%.1f", 100.0 * s.deposited / n) << "%, stopped_in_resist="
           << fmt("%.1f", 100.0 * s.in_mask / n) << "%\n";
    }
    return s;
  }

  const McImplantStats s = apply_mc_implant(st.mesh, silicon_mask(st), p, f);
  if (log) {
    const double n = static_cast<double>(p.ions);
    *log << "[implant] " << dop->symbol << " MC: E=" << fmt("%.4g", p.energy_kev)
         << " keV, deposited " << fmt("%.1f", 100.0 * s.deposited / n)
         << "% (Rp=" << fmt("%.4g", s.rp * 1e4) << " um)\n";
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

void mask(SimState& st, double x1, double x2, double y1, double y2,
          std::ostream* log) {
  need_mesh(st);
  if (!st.has_stack) throw std::runtime_error("no resist present; photo first");
  if (!(x2 > x1)) throw std::runtime_error("mask requires x2 > x1");
  if (!(y2 > y1)) throw std::runtime_error("mask requires y2 > y1");
  int opened = 0;
  const int nc_s = static_cast<int>(st.stack.cells.size());
  for (int ci = 0; ci < nc_s; ++ci) {
    if (st.stack_cell_mat[ci] != 1) continue;
    const Vec3& c = st.stack.cell_cent[ci];
    if (c.x >= x1 && c.x <= x2 && c.y >= y1 && c.y <= y2) {
      st.stack_cell_mat[ci] = 2;  // vacuum / developed opening
      ++opened;
    }
  }
  if (log)
    *log << "[mask] x=[" << fmt("%.4g", x1 * 1e4) << "," << fmt("%.4g", x2 * 1e4)
         << "] um -> opened " << opened << " cells\n";
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
    *log << "[diffuse] T=" << fmt("%.5g", opts.temp) << " K, time="
         << fmt("%.5g", opts.time) << " s, species:";
    for (const auto& fl : fields) *log << " " << fl.dopant->symbol;
    *log << "\n";
  }
  DiffusionSolver solver(st.mesh, silicon_mask(st), log);
  solver.run(fields, st.bcs, opts);
  st.last_temp = opts.temp;
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
    const double css = solid_solubility(*d, st.last_temp);
    extra.emplace_back(nc, 0.0);
    auto& act = extra.back();
    for (std::size_t i = 0; i < nc; ++i) {
      act[i] = (css > 0) ? std::min(conc[i], css) : conc[i];
      net[i] += (d->type == DopType::donor) ? act[i] : -act[i];
    }
  }
  std::size_t k = 0;
  for (const auto& [sym, conc] : st.fields) {
    (void)conc;
    scalars.emplace_back(sym + "_active", &extra[k++]);
  }
  scalars.emplace_back("NetDoping", &net);

  std::vector<int> region(st.mesh.cell_region.begin(), st.mesh.cell_region.end());
  std::vector<std::pair<std::string, const std::vector<int>*>> ints = {
      {"Region", &region}};
  write_vtu(path, st.mesh, scalars, ints);
  if (log) *log << "[save] wrote " << path << "\n";
}

}  // namespace proc
}  // namespace cp
