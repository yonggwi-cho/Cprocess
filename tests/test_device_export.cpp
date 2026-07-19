// Tests for P3-h: device-simulator export (VTU point data + meta.json).

#include <cmath>
#include <cstdio>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

#include "cprocess/deck.hpp"
#include "cprocess/device_export.hpp"
#include "cprocess/process.hpp"
#include "test_util.hpp"

using namespace cp;

static const char* kPrefix = "build_test_export";

static std::string slurp(const std::string& path) {
  std::ifstream in(path);
  std::ostringstream ss;
  ss << in.rdbuf();
  return ss.str();
}

static std::size_t count_substr(const std::string& hay, const std::string& needle) {
  std::size_t n = 0, pos = 0;
  while ((pos = hay.find(needle, pos)) != std::string::npos) {
    ++n;
    pos += needle.size();
  }
  return n;
}

// ---------------------------------------------------------------------------
// Test 1: proc::save regression -- no <PointData> in plain save() output.
// ---------------------------------------------------------------------------
static void test_save_no_point_data() {
  std::printf("test_save_no_point_data\n");
  std::ostringstream log;
  SimState st;
  proc::mesh_box(st, 0, 1.0e-4, 0, 1.0e-4, 0, 1.0e-4, 6, 6, 6, &log);
  proc::set_region(st, "silicon", -1, &log);
  proc::init(st, "B", 1e15, -1, &log);
  const std::string path = std::string(kPrefix) + "_save.vtu";
  proc::save(st, path, &log);
  std::string content = slurp(path);
  CHECK(content.find("<PointData>") == std::string::npos);
  std::remove(path.c_str());
}

static SimState make_flow_state() {
  std::ostringstream log;
  SimState st;
  proc::mesh_box(st, 0, 1.0e-4, 0, 1.0e-4, 0, 1.0e-4, 6, 6, 12, &log);
  proc::set_region(st, "silicon", -1, &log);
  proc::init(st, "B", 1e15, -1, &log);
  proc::implant_gauss(st, "As", 1e15, 30.0, 0, 0, 0, false, 0, 0, 0, 0, false,
                      "gauss", &log);
  DiffuseOpts opts;
  opts.temp = 1273.15;
  opts.time = 600.0;
  opts.verbosity = 0;
  proc::diffuse(st, opts, &log);
  return st;
}

// ---------------------------------------------------------------------------
// Test 2: well-formed PointData block.
// ---------------------------------------------------------------------------
static void test_point_data_wellformed() {
  std::printf("test_point_data_wellformed\n");
  std::ostringstream log;
  SimState st = make_flow_state();
  const std::size_t nn = st.mesh.nodes.size();

  proc::export_device(st, kPrefix, &log);
  std::string content = slurp(std::string(kPrefix) + ".vtu");

  std::size_t pd_open = count_substr(content, "<PointData>");
  std::size_t pd_close = count_substr(content, "</PointData>");
  CHECK(pd_open == 1);
  CHECK(pd_close == 1);
  std::size_t pos_pd = content.find("<PointData>");
  std::size_t pos_cd = content.find("<CellData>");
  CHECK(pos_pd != std::string::npos && pos_cd != std::string::npos);
  CHECK(pos_pd < pos_cd);

  for (const char* name : {"As", "B", "As_active", "B_active", "NetDoping"}) {
    std::string tag = std::string("Name=\"") + name + "\"";
    std::size_t p = content.find(tag);
    CHECK(p != std::string::npos && p < pos_cd);
  }

  std::size_t open_tags = count_substr(content, "<DataArray");
  std::size_t close_tags = count_substr(content, "</DataArray>");
  CHECK(open_tags == close_tags);

  // Count value lines in the PointData "As" array.
  std::string pd_block = content.substr(pos_pd, pos_cd - pos_pd);
  std::size_t as_tag = pd_block.find("Name=\"As\"");
  CHECK(as_tag != std::string::npos);
  std::size_t arr_start = pd_block.find('\n', as_tag) + 1;
  std::size_t arr_end = pd_block.find("</DataArray>", arr_start);
  std::string arr = pd_block.substr(arr_start, arr_end - arr_start);
  std::size_t nlines = count_substr(arr, "\n");
  CHECK(nlines == nn);

  std::remove((std::string(kPrefix) + ".vtu").c_str());
  std::remove((std::string(kPrefix) + ".meta.json").c_str());
}

// ---------------------------------------------------------------------------
// Test 3: averaging rule -- uniform field, and gas exclusion.
// ---------------------------------------------------------------------------
static void test_averaging_rule() {
  std::printf("test_averaging_rule\n");
  std::ostringstream log;
  SimState st;
  proc::mesh_box(st, 0, 1.0e-4, 0, 1.0e-4, 0, 1.0e-4, 4, 4, 4, &log);
  proc::set_region(st, "silicon", -1, &log);
  for (auto& v : st.fields["B"]) v = 3e17;
  if (st.fields.find("B") == st.fields.end())
    st.fields["B"] = std::vector<double>(st.mesh.cells.size(), 3e17);
  else
    st.fields["B"].assign(st.mesh.cells.size(), 3e17);

  std::vector<char> exclude(st.mesh.cells.size(), 0);
  std::vector<double> node_vals = cp::cell_to_node(st.mesh, st.fields["B"], exclude);
  for (double v : node_vals) CHECK_NEAR(v, 3e17, 3e17 * 1e-12);

  // gas exclusion: etch top half via polygon covering the whole footprint,
  // depth = half the box height -> top layer of cells becomes gas.
  SimState st2;
  proc::mesh_box(st2, 0, 1.0e-4, 0, 1.0e-4, 0, 1.0e-4, 4, 4, 4, &log);
  proc::set_region(st2, "silicon", -1, &log);
  st2.fields["B"].assign(st2.mesh.cells.size(), 5e17);
  std::vector<std::pair<double, double>> poly = {
      {-1e-4, -1e-4}, {2e-4, -1e-4}, {2e-4, 2e-4}, {-1e-4, 2e-4}};
  proc::etch(st2, 0.25e-4, poly, "", &log);

  std::vector<char> exclude2(st2.mesh.cells.size(), 0);
  for (std::size_t c = 0; c < st2.mesh.cells.size(); ++c) {
    auto it = st2.region_material.find(st2.mesh.cell_region[c]);
    std::string mat = (it != st2.region_material.end()) ? it->second : "silicon";
    exclude2[c] = (mat == "gas") ? 1 : 0;
  }
  bool any_excluded = false;
  for (char c : exclude2) if (c) any_excluded = true;
  CHECK(any_excluded);

  std::vector<double> node_vals2 = cp::cell_to_node(st2.mesh, st2.fields["B"], exclude2);
  // Node touching only gas cells -> 0. Find topmost node (max z) which
  // should only be adjacent to gas cells after the etch.
  double zmax = -1e300;
  for (const auto& p : st2.mesh.nodes) zmax = std::max(zmax, p.z);
  bool found_zero_top = false;
  for (std::size_t n = 0; n < st2.mesh.nodes.size(); ++n) {
    if (std::abs(st2.mesh.nodes[n].z - zmax) < 1e-12) {
      if (node_vals2[n] == 0.0) found_zero_top = true;
    }
  }
  CHECK(found_zero_top);
}

// ---------------------------------------------------------------------------
// Test 4: NetDoping sign convention + neutral exclusion + cell/save match.
// ---------------------------------------------------------------------------
static void test_net_doping_sign() {
  std::printf("test_net_doping_sign\n");
  std::ostringstream log;
  SimState st = make_flow_state();  // B background + As implant, diffused
  st.last_temp = 1273.15;

  proc::export_device(st, kPrefix, &log);

  // Compare cell NetDoping against proc::save's own construction.
  const std::size_t nc = st.mesh.cells.size();
  std::vector<double> net_ref(nc, 0.0);
  for (const auto& [sym, conc] : st.fields) {
    const Dopant* d = find_dopant(sym);
    if (!d) continue;
    for (std::size_t i = 0; i < nc; ++i) {
      double act = active_concentration(*d, conc[i], st.last_temp);
      if (d->type == DopType::neutral) continue;
      net_ref[i] += (d->type == DopType::donor) ? act : -act;
    }
  }

  std::string content = slurp(std::string(kPrefix) + ".vtu");
  // crude parse of CellData NetDoping array
  std::size_t cd_pos = content.find("<CellData>");
  std::string cd_block = content.substr(cd_pos);
  std::size_t nd_tag = cd_block.find("Name=\"NetDoping\"");
  CHECK(nd_tag != std::string::npos);
  std::size_t arr_start = cd_block.find('\n', nd_tag) + 1;
  std::size_t arr_end = cd_block.find("</DataArray>", arr_start);
  std::istringstream iss(cd_block.substr(arr_start, arr_end - arr_start));
  std::vector<double> net_vtu;
  double v;
  while (iss >> v) net_vtu.push_back(v);
  CHECK(net_vtu.size() == nc);
  for (std::size_t i = 0; i < nc; ++i) {
    double denom = std::max(std::abs(net_ref[i]), 1.0);
    CHECK_NEAR(net_vtu[i], net_ref[i], denom * 1e-6);  // %.6e ascii round-trip
  }

  // Point NetDoping sign checks.
  double zmax = -1e300;
  for (const auto& p : st.mesh.nodes) zmax = std::max(zmax, p.z);
  std::size_t pd_pos = content.find("<PointData>");
  std::size_t pd_end = content.find("<CellData>");
  std::string pd_block = content.substr(pd_pos, pd_end - pd_pos);
  std::size_t pnd_tag = pd_block.find("Name=\"NetDoping\"");
  CHECK(pnd_tag != std::string::npos);
  std::size_t parr_start = pd_block.find('\n', pnd_tag) + 1;
  std::size_t parr_end = pd_block.find("</DataArray>", parr_start);
  std::istringstream piss(pd_block.substr(parr_start, parr_end - parr_start));
  std::vector<double> net_point;
  while (piss >> v) net_point.push_back(v);
  CHECK(net_point.size() == st.mesh.nodes.size());

  bool near_surface_positive = false, deep_negative = false;
  for (std::size_t n = 0; n < st.mesh.nodes.size(); ++n) {
    double z = st.mesh.nodes[n].z;
    if (zmax - z < 0.1e-4) {
      if (net_point[n] > 0) near_surface_positive = true;
    }
    if (zmax - z > 0.5e-4) {
      CHECK_NEAR(net_point[n], -1e15, 1e15 * 0.5);
      deep_negative = true;
    }
  }
  CHECK(near_surface_positive);
  CHECK(deep_negative);

  std::remove((std::string(kPrefix) + ".vtu").c_str());
  std::remove((std::string(kPrefix) + ".meta.json").c_str());

  // Add a neutral species field and verify bit-identical NetDoping.
  st.fields["C"] = std::vector<double>(nc, 1e20);
  proc::export_device(st, kPrefix, &log);
  std::string content2 = slurp(std::string(kPrefix) + ".vtu");
  std::size_t cd2 = content2.find("<CellData>");
  std::string cd_block2 = content2.substr(cd2);
  std::size_t nd_tag2 = cd_block2.find("Name=\"NetDoping\"");
  std::size_t arr_start2 = cd_block2.find('\n', nd_tag2) + 1;
  std::size_t arr_end2 = cd_block2.find("</DataArray>", arr_start2);
  CHECK(cd_block2.substr(arr_start2, arr_end2 - arr_start2) ==
        cd_block.substr(arr_start, arr_end - arr_start));

  std::remove((std::string(kPrefix) + ".vtu").c_str());
  std::remove((std::string(kPrefix) + ".meta.json").c_str());
}

// ---------------------------------------------------------------------------
// Test 5: meta.json schema conformance.
// ---------------------------------------------------------------------------
static void test_meta_json() {
  std::printf("test_meta_json\n");
  std::ostringstream log;
  SimState st = make_flow_state();
  proc::deposit(st, "oxide", 0.1e-4, 2, {}, &log);
  st.fields["C"] = std::vector<double>(st.mesh.cells.size(), 1e18);

  proc::export_device(st, kPrefix, &log);
  std::string meta = slurp(std::string(kPrefix) + ".meta.json");

  CHECK(meta.find("\"format\": \"cprocess-device\"") != std::string::npos);
  CHECK(meta.find("\"version\": 1") != std::string::npos);
  CHECK(meta.find("\"convention\": \"ND-NA\"") != std::string::npos);
  CHECK(meta.find("\"excluded_types\": [\"neutral\"]") != std::string::npos);
  CHECK(meta.find("\"material\": \"oxide\"") != std::string::npos);
  std::size_t c_pos = meta.find("\"symbol\": \"C\",");
  CHECK(c_pos != std::string::npos);
  std::size_t line_end = meta.find('\n', c_pos);
  std::size_t line_start = meta.rfind('\n', c_pos);
  std::string line = meta.substr(line_start, line_end - line_start);
  CHECK(line.find("\"in_net_doping\": false") != std::string::npos);
  for (const char* p : {"xmin", "xmax", "ymin", "ymax", "zmin", "zmax"})
    CHECK(meta.find(p) != std::string::npos);
  CHECK(meta.find("\"length\": \"um\"") != std::string::npos);

  std::remove((std::string(kPrefix) + ".vtu").c_str());
  std::remove((std::string(kPrefix) + ".meta.json").c_str());
}

// ---------------------------------------------------------------------------
// Test 6: error cases.
// ---------------------------------------------------------------------------
static void test_errors() {
  std::printf("test_errors\n");
  std::ostringstream log;
  SimState st;
  proc::mesh_box(st, 0, 1.0e-4, 0, 1.0e-4, 0, 1.0e-4, 4, 4, 8, &log);
  proc::set_region(st, "silicon", -1, &log);
  proc::init(st, "B", 1e15, -1, &log);

  proc::photo(st, 0.3e-4, 2, &log);
  // W-7: a live resist stack no longer throws — warn and export without it.
  std::ostringstream wlog;
  bool threw = false;
  try {
    proc::export_device(st, kPrefix, &wlog);
  } catch (const std::exception&) {
    threw = true;
  }
  CHECK(!threw);
  CHECK(wlog.str().find("warning") != std::string::npos);

  proc::strip(st, &log);
  threw = false;
  try {
    proc::export_device(st, kPrefix, &log);
  } catch (const std::exception&) {
    threw = true;
  }
  CHECK(!threw);
  std::remove((std::string(kPrefix) + ".vtu").c_str());
  std::remove((std::string(kPrefix) + ".meta.json").c_str());

  threw = false;
  try {
    proc::export_device(st, "/nonexistent_dir/x", &log);
  } catch (const std::exception&) {
    threw = true;
  }
  CHECK(threw);

  SimState empty;
  threw = false;
  try {
    proc::export_device(empty, kPrefix, &log);
  } catch (const std::exception&) {
    threw = true;
  }
  CHECK(threw);
}

int main() {
  test_save_no_point_data();
  test_point_data_wellformed();
  test_averaging_rule();
  test_net_doping_sign();
  test_meta_json();
  test_errors();
  std::printf("test_device_export: all tests passed\n");
  return 0;
}
