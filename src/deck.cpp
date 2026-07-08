#include "cprocess/deck.hpp"

#include <algorithm>
#include <cctype>
#include <cstdlib>
#include <istream>
#include <map>
#include <ostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include "cprocess/process.hpp"

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

// A mask window must be given as all of x1,x2,y1,y2 — or none.
bool parse_window(const Cmd& c, double& x1, double& x2, double& y1, double& y2) {
  const int nw = c.has("x1") + c.has("x2") + c.has("y1") + c.has("y2");
  if (nw == 0) return false;
  if (nw != 4) c.fail("mask window needs all of x1=, x2=, y1=, y2=");
  x1 = c.num("x1", Unit::length);
  x2 = c.num("x2", Unit::length);
  y1 = c.num("y1", Unit::length);
  y2 = c.num("y2", Unit::length);
  if (!(x2 > x1) || !(y2 > y1)) c.fail("window must have x2>x1, y2>y1");
  return true;
}

void cmd_mesh(SimState& st, const Cmd& c, std::ostream& log) {
  if (c.bare.size() < 2) c.fail("expected 'mesh box ...' or 'mesh gmsh ...'");
  const std::string sub = lower(c.bare[1]);
  if (sub == "box") {
    proc::mesh_box(st,
        c.num_or("xmin", Unit::length, 0), c.num("xmax", Unit::length),
        c.num_or("ymin", Unit::length, 0), c.num("ymax", Unit::length),
        c.num_or("zmin", Unit::length, 0), c.num("zmax", Unit::length),
        static_cast<int>(c.num("nx", Unit::none)),
        static_cast<int>(c.num("ny", Unit::none)),
        static_cast<int>(c.num("nz", Unit::none)), &log);
  } else if (sub == "gmsh") {
    proc::mesh_gmsh(st, c.str("file"), c.num_or("scale", Unit::length, 1.0), &log);
  } else {
    c.fail("unknown mesh type '" + sub + "'");
  }
}

void cmd_region(SimState& st, const Cmd& c, std::ostream& log) {
  const std::string mat = c.str("material");
  if (!c.bare.empty() && c.bare.size() >= 2 && lower(c.bare[1]) == "all") {
    proc::set_region(st, mat, -1, &log);
  } else if (c.has("tag")) {
    proc::set_region(st, mat, static_cast<int>(c.num("tag", Unit::none)), &log);
  } else if (c.has("name")) {
    proc::set_region(st, mat, proc::resolve_region(st, c.str("name")), &log);
  } else {
    c.fail("expected tag=, name= or 'region all'");
  }
}

void cmd_init(SimState& st, const Cmd& c, std::ostream& log) {
  const int region = c.has("region")
      ? proc::resolve_region(st, c.str("region")) : -1;
  proc::init(st, c.str("species"), c.num("conc", Unit::none), region, &log);
}

void cmd_implant(SimState& st, const Cmd& c, std::ostream& log) {
  const std::string species = c.str("species");
  const double dose = c.num("dose", Unit::none);
  double x1 = 0, x2 = 0, y1 = 0, y2 = 0;
  const bool has_window = parse_window(c, x1, x2, y1, y2);
  const std::string method =
      c.has("method") ? lower(c.kv.at("method")) : std::string("gauss");

  if (method == "mc" || method == "montecarlo") {
    proc::implant_mc(st, species, dose, c.num("energy", Unit::energy),
        static_cast<long long>(c.num_or("ions", Unit::none, 100000)),
        c.num_or("tilt", Unit::none, 0), c.num_or("rotation", Unit::none, 0),
        static_cast<unsigned long long>(c.num_or("seed", Unit::none, 1)),
        static_cast<int>(c.num_or("threads", Unit::none, 0)),
        c.flag_or("channeling", false), has_window, x1, x2, y1, y2,
        c.flag_or("lateral_wrap", false), c.flag_or("damage", false), &log);
    return;
  }
  if (method != "gauss" && method != "gaussian" && method != "analytic")
    c.fail("unknown method '" + method + "' (gauss or mc)");

  double rp = 0, drp = 0, energy = 0;
  if (c.has("rp") || c.has("drp")) {
    rp = c.num("rp", Unit::length);
    drp = c.num("drp", Unit::length);
  } else {
    energy = c.num("energy", Unit::energy);
  }
  const std::string profile =
      c.has("profile") ? lower(c.kv.at("profile")) : std::string("gauss");
  proc::implant_gauss(st, species, dose, energy, rp, drp,
                      c.num_or("drl", Unit::length, 0), has_window, x1, x2, y1, y2,
                      c.flag_or("damage", false), profile, &log);
}

void cmd_photo(SimState& st, const Cmd& c, std::ostream& log) {
  proc::photo(st, c.num("resist", Unit::length),
              static_cast<int>(c.num_or("nz", Unit::none, 4)), &log);
}

void cmd_mask(SimState& st, const Cmd& c, std::ostream& log) {
  const BBox bb = st.has_stack ? st.stack.bbox() : st.mesh.bbox();
  proc::mask(st, c.num("x1", Unit::length), c.num("x2", Unit::length),
             c.num_or("y1", Unit::length, bb.lo.y),
             c.num_or("y2", Unit::length, bb.hi.y), &log);
}

void cmd_strip(SimState& st, const Cmd&, std::ostream& log) {
  proc::strip(st, &log);
}

void cmd_bc(SimState& st, const Cmd& c, std::ostream& log) {
  if (c.bare.size() >= 2 && lower(c.bare[1]) == "clear") {
    proc::clear_bc(st, &log);
    return;
  }
  const int patch = st.mesh.find_patch(c.str("patch"));
  if (patch < 0) c.fail("unknown patch '" + c.str("patch") + "'");
  proc::add_bc(st, c.str("species"), patch, c.num("conc", Unit::none), &log);
}

void cmd_diffuse(SimState& st, const Cmd& c, std::ostream& log) {
  DiffuseOpts o;
  o.time = c.num("time", Unit::time);
  o.temp = c.num("temp", Unit::temp);
  o.dt = c.num_or("dt", Unit::time, 0);
  o.field_enh = c.flag_or("fieldenh", true);
  o.nonortho = c.flag_or("nonortho", true);
  if (o.time <= 0) c.fail("time must be > 0");
  if (c.flag_or("ted", false))
    proc::diffuse_ted(st, o, &log);
  else
    proc::diffuse(st, o, &log);
}

void cmd_oxidize(SimState& st, const Cmd& c, std::ostream& log) {
  const double time_s = c.num("time", Unit::time);
  const double temp_k = c.num("temp", Unit::temp);
  std::string ambient = c.has("ambient") ? lower(c.str("ambient")) : "dry";
  if (ambient != "dry" && ambient != "wet")
    c.fail("ambient must be 'dry' or 'wet'");
  proc::oxidize(st, time_s, temp_k, ambient == "wet", &log);
}

void cmd_save(SimState& st, const Cmd& c, std::ostream& log) {
  proc::save(st, c.has("file") ? c.str("file") : "out.vtu", &log);
}

void cmd_print(SimState& st, const Cmd& c, std::ostream& log) {
  if (!st.has_mesh) c.fail("no mesh defined yet");
  const BBox b = st.mesh.bbox();
  log << "[print] mesh: " << st.mesh.cells.size() << " tets\n";
  for (const auto& [sym, conc] : st.fields) {
    double mass = 0, peak = 0;
    for (std::size_t i = 0; i < conc.size(); ++i) {
      mass += conc[i] * st.mesh.cell_vol[i];
      if (conc[i] > peak) peak = conc[i];
    }
    const double area = (b.hi.x - b.lo.x) * (b.hi.y - b.lo.y);
    char buf[256];
    std::snprintf(buf, sizeof(buf),
        "[print] %s: integral=%.4g atoms (%.4g cm^-2 avg), peak=%.4g cm^-3\n",
        sym.c_str(), mass, mass / area, peak);
    log << buf;
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
    else if (c.name == "photo") cmd_photo(st, c, log);
    else if (c.name == "mask") cmd_mask(st, c, log);
    else if (c.name == "strip") cmd_strip(st, c, log);
    else if (c.name == "bc") cmd_bc(st, c, log);
    else if (c.name == "diffuse") cmd_diffuse(st, c, log);
    else if (c.name == "oxidize") cmd_oxidize(st, c, log);
    else if (c.name == "save") cmd_save(st, c, log);
    else if (c.name == "print") cmd_print(st, c, log);
    else if (c.name == "stop") break;
    else c.fail("unknown command");
  }
}

}  // namespace cp
