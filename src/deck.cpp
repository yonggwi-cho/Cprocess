#include "cprocess/deck.hpp"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <istream>
#include <ostream>
#include <sstream>
#include <stdexcept>

#include "cprocess/gmsh_reader.hpp"
#include "cprocess/implant.hpp"
#include "cprocess/vtk_writer.hpp"

namespace cp {

namespace {

enum class Unit { none, length, time, temp, energy };

std::string lower(std::string s) {
  std::transform(s.begin(), s.end(), s.begin(),
                 [](unsigned char c) { return std::tolower(c); });
  return s;
}

// Parses "0.5um", "30min", "1000C", "80keV", "1e15". Defaults (no suffix):
// cm, s, K, keV.
bool parse_quantity(const std::string& text, Unit unit, double& out) {
  const char* cs = text.c_str();
  char* end = nullptr;
  const double v = std::strtod(cs, &end);
  if (end == cs) return false;
  std::string suf = lower(std::string(end));
  suf.erase(std::remove_if(suf.begin(), suf.end(),
                           [](unsigned char c) { return std::isspace(c); }),
            suf.end());
  if (suf.empty()) {
    out = v;
    return true;
  }
  switch (unit) {
    case Unit::length:
      if (suf == "cm") out = v;
      else if (suf == "mm") out = v * 0.1;
      else if (suf == "um" || suf == "\xc2\xb5m") out = v * 1e-4;
      else if (suf == "nm") out = v * 1e-7;
      else if (suf == "m") out = v * 100.0;
      else return false;
      return true;
    case Unit::time:
      if (suf == "s" || suf == "sec") out = v;
      else if (suf == "ms") out = v * 1e-3;
      else if (suf == "min") out = v * 60.0;
      else if (suf == "h" || suf == "hr") out = v * 3600.0;
      else return false;
      return true;
    case Unit::temp:
      if (suf == "k") out = v;
      else if (suf == "c") out = v + 273.15;
      else return false;
      return true;
    case Unit::energy:
      if (suf == "kev") out = v;
      else if (suf == "ev") out = v * 1e-3;
      else if (suf == "mev") out = v * 1e3;
      else return false;
      return true;
    case Unit::none:
      return false;
  }
  return false;
}

struct Cmd {
  std::string name;
  std::vector<std::string> bare;
  std::map<std::string, std::string> kv;
  int line = 0;

  [[noreturn]] void fail(const std::string& msg) const {
    throw std::runtime_error("deck line " + std::to_string(line) + " (" + name +
                             "): " + msg);
  }
  bool has(const std::string& k) const { return kv.count(k) > 0; }
  std::string str(const std::string& k) const {
    auto it = kv.find(k);
    if (it == kv.end()) fail("missing argument '" + k + "='");
    return it->second;
  }
  double num(const std::string& k, Unit u) const {
    double v;
    if (!parse_quantity(str(k), u, v)) fail("cannot parse '" + k + "=" + str(k) + "'");
    return v;
  }
  double num_or(const std::string& k, Unit u, double dflt) const {
    if (!has(k)) return dflt;
    return num(k, u);
  }
  bool flag_or(const std::string& k, bool dflt) const {
    if (!has(k)) return dflt;
    const std::string v = lower(kv.at(k));
    if (v == "on" || v == "true" || v == "1" || v == "yes") return true;
    if (v == "off" || v == "false" || v == "0" || v == "no") return false;
    fail("expected on/off for '" + k + "='");
  }
};

bool is_silicon(const std::string& mat) {
  const std::string m = lower(mat);
  return m == "silicon" || m == "si";
}

std::vector<char> silicon_mask(const SimState& st) {
  std::vector<char> mask(st.mesh.cells.size(), 1);
  for (std::size_t i = 0; i < st.mesh.cells.size(); ++i) {
    auto it = st.region_material.find(st.mesh.cell_region[i]);
    mask[i] = (it == st.region_material.end()) || is_silicon(it->second);
  }
  return mask;
}

int region_by_name_or_tag(const SimState& st, const Cmd& c, const std::string& v) {
  char* end = nullptr;
  const long tag = std::strtol(v.c_str(), &end, 10);
  if (end != v.c_str() && *end == '\0') return static_cast<int>(tag);
  for (const auto& [t, n] : st.mesh.region_names)
    if (lower(n) == lower(v)) return t;
  c.fail("unknown region '" + v + "'");
}

const Dopant* dopant_arg(const Cmd& c) {
  const Dopant* d = find_dopant(c.str("species"));
  if (!d) c.fail("unknown species '" + c.str("species") + "' (B, P, As, Sb)");
  return d;
}

void need_mesh(const SimState& st, const Cmd& c) {
  if (!st.has_mesh) c.fail("no mesh defined yet (use 'mesh box' or 'mesh gmsh')");
}

std::string fmt(const char* f, double v) {
  char buf[64];
  std::snprintf(buf, sizeof(buf), f, v);
  return buf;
}

void cmd_mesh(SimState& st, const Cmd& c, std::ostream& log) {
  if (c.bare.size() < 2) c.fail("expected 'mesh box ...' or 'mesh gmsh ...'");
  const std::string sub = lower(c.bare[1]);
  if (sub == "box") {
    st.mesh = make_box_mesh(
        c.num_or("xmin", Unit::length, 0), c.num("xmax", Unit::length),
        c.num_or("ymin", Unit::length, 0), c.num("ymax", Unit::length),
        c.num_or("zmin", Unit::length, 0), c.num("zmax", Unit::length),
        static_cast<int>(c.num("nx", Unit::none)),
        static_cast<int>(c.num("ny", Unit::none)),
        static_cast<int>(c.num("nz", Unit::none)));
  } else if (sub == "gmsh") {
    st.mesh = read_gmsh(c.str("file"), c.num_or("scale", Unit::length, 1.0), &log);
  } else {
    c.fail("unknown mesh type '" + sub + "'");
  }
  st.has_mesh = true;
  st.fields.clear();
  st.bcs.clear();
  st.region_material.clear();
  for (int tag : st.mesh.region_tags()) st.region_material[tag] = "silicon";

  const BBox b = st.mesh.bbox();
  log << "[mesh] " << st.mesh.cells.size() << " tets, " << st.mesh.nodes.size()
      << " nodes, extent " << fmt("%.3g", (b.hi.x - b.lo.x) * 1e4) << " x "
      << fmt("%.3g", (b.hi.y - b.lo.y) * 1e4) << " x "
      << fmt("%.3g", (b.hi.z - b.lo.z) * 1e4) << " um, min orthogonality "
      << fmt("%.3f", st.mesh.min_orthogonality()) << "\n[mesh] patches:";
  for (const auto& p : st.mesh.patch_names) log << " " << p;
  log << "\n[mesh] regions:";
  for (int tag : st.mesh.region_tags()) {
    log << " " << tag;
    auto it = st.mesh.region_names.find(tag);
    if (it != st.mesh.region_names.end() && !it->second.empty())
      log << "(" << it->second << ")";
  }
  log << "\n";
}

void cmd_region(SimState& st, const Cmd& c, std::ostream& log) {
  need_mesh(st, c);
  const std::string mat = lower(c.str("material"));
  static const char* known[] = {"silicon", "si", "oxide", "nitride", "poly",
                                "polysilicon", "gas"};
  if (std::none_of(std::begin(known), std::end(known),
                   [&](const char* k) { return mat == k; }))
    c.fail("unknown material '" + mat + "'");
  std::vector<int> tags;
  if (!c.bare.empty() && c.bare.size() >= 2 && lower(c.bare[1]) == "all") {
    tags = st.mesh.region_tags();
  } else if (c.has("tag")) {
    tags.push_back(static_cast<int>(c.num("tag", Unit::none)));
  } else if (c.has("name")) {
    tags.push_back(region_by_name_or_tag(st, c, c.str("name")));
  } else {
    c.fail("expected tag=, name= or 'region all'");
  }
  for (int t : tags) st.region_material[t] = mat;
  log << "[region]";
  for (int t : tags) log << " " << t;
  log << " -> " << mat << "\n";
}

void cmd_init(SimState& st, const Cmd& c, std::ostream& log) {
  need_mesh(st, c);
  const Dopant* d = dopant_arg(c);
  const double conc = c.num("conc", Unit::none);
  if (conc < 0) c.fail("conc must be >= 0");
  int region = -1;
  if (c.has("region")) region = region_by_name_or_tag(st, c, c.str("region"));
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
  log << "[init] " << d->symbol << " = " << fmt("%.3g", conc) << " cm^-3 in "
      << nset << " cells\n";
}

void cmd_implant(SimState& st, const Cmd& c, std::ostream& log) {
  need_mesh(st, c);
  ImplantParams p;
  p.dopant = dopant_arg(c);
  p.dose = c.num("dose", Unit::none);
  if (c.has("rp") || c.has("drp")) {
    p.rp = c.num("rp", Unit::length);
    p.drp = c.num("drp", Unit::length);
  } else {
    const double e = c.num("energy", Unit::energy);
    if (!implant_range(*p.dopant, e, p.rp, p.drp))
      c.fail("no range table for species; give rp= and drp=");
    const auto& tab = p.dopant->range;
    if (e < tab.front()[0] || e > tab.back()[0])
      log << "[implant] warning: energy outside range table ("
          << tab.front()[0] << "-" << tab.back()[0] << " keV), clamped\n";
  }
  p.drl = c.num_or("drl", Unit::length, 0);
  const int nw = c.has("x1") + c.has("x2") + c.has("y1") + c.has("y2");
  if (nw == 4) {
    p.has_window = true;
    p.x1 = c.num("x1", Unit::length);
    p.x2 = c.num("x2", Unit::length);
    p.y1 = c.num("y1", Unit::length);
    p.y2 = c.num("y2", Unit::length);
    if (!(p.x2 > p.x1) || !(p.y2 > p.y1)) c.fail("window must have x2>x1, y2>y1");
  } else if (nw != 0) {
    c.fail("mask window needs all of x1=, x2=, y1=, y2=");
  }

  auto& f = st.fields[p.dopant->symbol];
  f.resize(st.mesh.cells.size(), 0.0);
  const double atoms = apply_implant(st.mesh, silicon_mask(st), p, f);
  const BBox b = st.mesh.bbox();
  const double area = p.has_window ? (p.x2 - p.x1) * (p.y2 - p.y1)
                                   : (b.hi.x - b.lo.x) * (b.hi.y - b.lo.y);
  log << "[implant] " << p.dopant->symbol << " dose=" << fmt("%.3g", p.dose)
      << " cm^-2, Rp=" << fmt("%.4g", p.rp * 1e4)
      << " um, dRp=" << fmt("%.4g", p.drp * 1e4) << " um";
  if (p.has_window)
    log << ", window [" << fmt("%.3g", p.x1 * 1e4) << ","
        << fmt("%.3g", p.x2 * 1e4) << "]x[" << fmt("%.3g", p.y1 * 1e4) << ","
        << fmt("%.3g", p.y2 * 1e4) << "] um";
  log << " -> integrated " << fmt("%.4g", atoms / area) << " cm^-2\n";
}

void cmd_bc(SimState& st, const Cmd& c, std::ostream& log) {
  need_mesh(st, c);
  if (c.bare.size() >= 2 && lower(c.bare[1]) == "clear") {
    st.bcs.clear();
    log << "[bc] cleared\n";
    return;
  }
  DirichletBC bc;
  const Dopant* d = dopant_arg(c);
  bc.species = d->symbol;
  bc.patch = st.mesh.find_patch(c.str("patch"));
  if (bc.patch < 0) c.fail("unknown patch '" + c.str("patch") + "'");
  bc.conc = c.num("conc", Unit::none);
  if (bc.conc < 0) c.fail("conc must be >= 0");
  st.bcs.push_back(bc);
  st.fields[d->symbol].resize(st.mesh.cells.size(), 0.0);
  log << "[bc] " << bc.species << " = " << fmt("%.3g", bc.conc) << " cm^-3 on '"
      << c.str("patch") << "'\n";
}

void cmd_diffuse(SimState& st, const Cmd& c, std::ostream& log) {
  need_mesh(st, c);
  DiffuseOpts o;
  o.time = c.num("time", Unit::time);
  o.temp = c.num("temp", Unit::temp);
  o.dt = c.num_or("dt", Unit::time, 0);
  o.field_enh = c.flag_or("fieldenh", true);
  o.nonortho = c.flag_or("nonortho", true);
  if (o.time <= 0) c.fail("time must be > 0");
  if (o.temp < 600 || o.temp > 1800)
    log << "[diffuse] warning: T=" << o.temp
        << " K outside the usual 600-1800 K model range\n";

  std::vector<SpeciesField> fields;
  for (auto& [sym, conc] : st.fields) {
    const Dopant* d = find_dopant(sym);
    if (d) fields.push_back({d, &conc});
  }
  if (fields.empty()) {
    log << "[diffuse] no dopants present, nothing to do\n";
    return;
  }
  log << "[diffuse] T=" << fmt("%.5g", o.temp) << " K, time=" << fmt("%.5g", o.time)
      << " s, dt=" << fmt("%.4g", o.dt > 0 ? std::min(o.dt, o.time) : o.time / 50)
      << " s, ni=" << fmt("%.3g", ni_si(o.temp)) << " cm^-3, species:";
  for (const auto& f : fields) log << " " << f.dopant->symbol;
  log << "\n";

  DiffusionSolver solver(st.mesh, silicon_mask(st), &log);
  solver.run(fields, st.bcs, o);
  st.last_temp = o.temp;
}

void cmd_save(SimState& st, const Cmd& c, std::ostream& log) {
  need_mesh(st, c);
  const std::string path = c.has("file") ? c.str("file") : "out.vtu";
  const std::size_t nc = st.mesh.cells.size();

  std::vector<std::vector<double>> extra;
  std::vector<std::pair<std::string, const std::vector<double>*>> scalars;
  for (const auto& [sym, conc] : st.fields) scalars.emplace_back(sym, &conc);

  // Electrically active concentrations: clamped at the solid solubility of
  // the last anneal temperature.
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
  log << "[save] wrote " << path << " (" << scalars.size() << " scalar fields, "
      << "active clamp at T=" << fmt("%.5g", st.last_temp) << " K)\n";
}

void cmd_print(SimState& st, const Cmd& c, std::ostream& log) {
  need_mesh(st, c);
  const BBox b = st.mesh.bbox();
  log << "[print] mesh: " << st.mesh.cells.size() << " tets, volume "
      << fmt("%.6g", st.mesh.total_volume()) << " cm^3\n";
  for (const auto& [sym, conc] : st.fields) {
    double mass = 0, peak = 0;
    Vec3 ploc{};
    for (std::size_t i = 0; i < conc.size(); ++i) {
      mass += conc[i] * st.mesh.cell_vol[i];
      if (conc[i] > peak) {
        peak = conc[i];
        ploc = st.mesh.cell_cent[i];
      }
    }
    const double area = (b.hi.x - b.lo.x) * (b.hi.y - b.lo.y);
    log << "[print] " << sym << ": integral=" << fmt("%.4g", mass)
        << " atoms (" << fmt("%.4g", mass / area) << " cm^-2 avg), peak="
        << fmt("%.4g", peak) << " cm^-3 at (" << fmt("%.3g", ploc.x * 1e4)
        << ", " << fmt("%.3g", ploc.y * 1e4) << ", " << fmt("%.3g", ploc.z * 1e4)
        << ") um\n";
  }
}

}  // namespace

void run_deck(std::istream& in, SimState& st, std::ostream& log) {
  std::string line;
  int lineno = 0;
  while (std::getline(in, line)) {
    ++lineno;
    const auto hash = line.find('#');
    if (hash != std::string::npos) line.erase(hash);

    Cmd c;
    c.line = lineno;
    std::istringstream is(line);
    std::string tok;
    while (is >> tok) {
      const auto eq = tok.find('=');
      if (c.bare.empty()) {
        c.bare.push_back(tok);
      } else if (eq != std::string::npos && eq > 0) {
        c.kv[lower(tok.substr(0, eq))] = tok.substr(eq + 1);
      } else {
        c.bare.push_back(tok);
      }
    }
    if (c.bare.empty()) continue;
    c.name = lower(c.bare[0]);

    if (c.name == "mesh") cmd_mesh(st, c, log);
    else if (c.name == "region") cmd_region(st, c, log);
    else if (c.name == "init") cmd_init(st, c, log);
    else if (c.name == "implant") cmd_implant(st, c, log);
    else if (c.name == "bc") cmd_bc(st, c, log);
    else if (c.name == "diffuse") cmd_diffuse(st, c, log);
    else if (c.name == "save") cmd_save(st, c, log);
    else if (c.name == "print") cmd_print(st, c, log);
    else if (c.name == "stop") break;
    else c.fail("unknown command");
  }
}

}  // namespace cp
