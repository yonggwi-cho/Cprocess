#include "cprocess/device_export.hpp"

#include <cstdio>
#include <fstream>
#include <map>
#include <set>
#include <sstream>
#include <stdexcept>

#include "cprocess/materials.hpp"
#include "cprocess/vtk_writer.hpp"

namespace cp {

std::vector<double> cell_to_node(const Mesh& mesh,
                                  const std::vector<double>& f,
                                  const std::vector<char>& exclude) {
  const std::size_t nn = mesh.nodes.size();
  std::vector<double> num(nn, 0.0), den(nn, 0.0);
  for (std::size_t c = 0; c < mesh.cells.size(); ++c) {
    if (c < exclude.size() && exclude[c]) continue;
    const double v = c < mesh.cell_vol.size() ? mesh.cell_vol[c] : 0.0;
    const double val = c < f.size() ? f[c] : 0.0;
    for (int n : mesh.cells[c]) {
      num[static_cast<std::size_t>(n)] += v * val;
      den[static_cast<std::size_t>(n)] += v;
    }
  }
  std::vector<double> out(nn, 0.0);
  for (std::size_t n = 0; n < nn; ++n)
    if (den[n] > 0.0) out[n] = num[n] / den[n];
  return out;
}

namespace {

std::string json_escape(const std::string& s) {
  std::string out;
  out.reserve(s.size());
  for (char c : s) {
    if (c == '"' || c == '\\') out.push_back('\\');
    out.push_back(c);
  }
  return out;
}

std::string dop_type_name(DopType t) {
  switch (t) {
    case DopType::donor: return "donor";
    case DopType::acceptor: return "acceptor";
    case DopType::neutral: return "neutral";
  }
  return "neutral";
}

bool is_gas_cell(const SimState& st, int region_tag) {
  auto it = st.region_material.find(region_tag);
  const std::string mat = (it != st.region_material.end()) ? it->second : "silicon";
  return material_id(mat) == kMatGas;
}

}  // namespace

void write_device(const SimState& st, const std::string& prefix) {
  const std::size_t nc = st.mesh.cells.size();
  const std::size_t nn = st.mesh.nodes.size();

  // ── CellData: identical construction to proc::save ──
  std::vector<std::pair<std::string, const std::vector<double>*>> scalars;
  for (const auto& [sym, conc] : st.fields) scalars.emplace_back(sym, &conc);

  std::vector<std::vector<double>> extra;
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
  std::vector<std::string> cell_field_names;
  for (const auto& [sym, conc] : st.fields) {
    (void)conc;
    cell_field_names.push_back(sym);
  }
  for (const auto& [sym, conc] : st.fields) {
    (void)conc;
    if (!find_dopant(sym)) continue;
    scalars.emplace_back(sym + "_active", &extra[k++]);
    cell_field_names.push_back(sym + "_active");
  }
  scalars.emplace_back("NetDoping", &net);
  cell_field_names.push_back("NetDoping");

  std::vector<int> region(st.mesh.cell_region.begin(), st.mesh.cell_region.end());
  std::vector<std::pair<std::string, const std::vector<int>*>> ints = {
      {"Region", &region}};
  cell_field_names.push_back("Region");

  // ── PointData: dopant fields, their _active, and NetDoping, node-averaged ──
  std::vector<char> exclude(nc, 0);
  for (std::size_t c = 0; c < nc; ++c)
    exclude[c] = is_gas_cell(st, st.mesh.cell_region[c]) ? 1 : 0;

  std::vector<std::vector<double>> point_data;
  std::vector<std::pair<std::string, const std::vector<double>*>> point_scalars;
  std::vector<std::string> point_field_names;

  // Reserve capacity so pointers into point_data stay valid.
  std::size_t n_point_arrays = 1;  // NetDoping
  for (const auto& [sym, conc] : st.fields) {
    (void)conc;
    if (find_dopant(sym)) n_point_arrays += 2;  // field + _active
  }
  point_data.reserve(n_point_arrays);

  k = 0;
  for (const auto& [sym, conc] : st.fields) {
    const Dopant* d = find_dopant(sym);
    if (!d) continue;
    point_data.push_back(cell_to_node(st.mesh, conc, exclude));
    point_scalars.emplace_back(sym, &point_data.back());
    point_field_names.push_back(sym);
    point_data.push_back(cell_to_node(st.mesh, extra[k++], exclude));
    point_scalars.emplace_back(sym + "_active", &point_data.back());
    point_field_names.push_back(sym + "_active");
  }
  point_data.push_back(cell_to_node(st.mesh, net, exclude));
  point_scalars.emplace_back("NetDoping", &point_data.back());
  point_field_names.push_back("NetDoping");

  const std::string vtu_path = prefix + ".vtu";
  write_vtu(vtu_path, st.mesh, scalars, ints, point_scalars);

  // ── meta.json ──
  const std::string meta_path = prefix + ".meta.json";
  std::ofstream out(meta_path);
  if (!out) throw std::runtime_error("write_device: cannot open for writing: " + meta_path);

  // basename of vtu_path for "files.vtu"
  std::string vtu_base = vtu_path;
  auto slash = vtu_base.find_last_of("/\\");
  if (slash != std::string::npos) vtu_base = vtu_base.substr(slash + 1);

  char buf[64];

  out << "{\n";
  out << "  \"format\": \"cprocess-device\",\n";
  out << "  \"version\": 1,\n";
  out << "  \"generator\": \"cprocess proc::export_device (P3-h)\",\n";
  out << "  \"files\": {\n";
  out << "    \"vtu\": \"" << json_escape(vtu_base) << "\"\n";
  out << "  },\n";
  out << "  \"units\": {\n";
  out << "    \"length\": \"um\",\n";
  out << "    \"concentration\": \"cm^-3\",\n";
  out << "    \"temperature\": \"K\",\n";
  out << "    \"note\": \"VTU point coordinates are in micrometres (write_vtu "
         "to_um=1e4); all concentration arrays are cm^-3\"\n";
  out << "  },\n";
  out << "  \"mesh\": {\n";
  out << "    \"n_nodes\": " << nn << ",\n";
  out << "    \"n_cells\": " << nc << ",\n";
  out << "    \"cell_type\": \"tetrahedron\"\n";
  out << "  },\n";

  // regions: union of region_names and region_material tags, sorted.
  std::set<int> tags;
  for (const auto& [tag, name] : st.mesh.region_names) { (void)name; tags.insert(tag); }
  for (const auto& [tag, mat] : st.region_material) { (void)mat; tags.insert(tag); }
  out << "  \"regions\": [\n";
  {
    std::size_t i = 0, ntags = tags.size();
    for (int tag : tags) {
      auto nit = st.mesh.region_names.find(tag);
      std::string name = (nit != st.mesh.region_names.end()) ? nit->second
                                                              : ("region" + std::to_string(tag));
      auto mit = st.region_material.find(tag);
      std::string mat = (mit != st.region_material.end()) ? mit->second : "silicon";
      out << "    { \"tag\": " << tag << ", \"name\": \"" << json_escape(name)
          << "\", \"material\": \"" << json_escape(mat) << "\" }";
      out << (++i < ntags ? ",\n" : "\n");
    }
  }
  out << "  ],\n";

  out << "  \"boundaries\": [";
  for (std::size_t i = 0; i < st.mesh.patch_names.size(); ++i) {
    out << "\"" << json_escape(st.mesh.patch_names[i]) << "\"";
    if (i + 1 < st.mesh.patch_names.size()) out << ", ";
  }
  out << "],\n";

  out << "  \"species\": [\n";
  {
    std::vector<std::string> syms;
    for (const auto& [sym, conc] : st.fields) {
      (void)conc;
      if (find_dopant(sym)) syms.push_back(sym);
    }
    for (std::size_t i = 0; i < syms.size(); ++i) {
      const Dopant* d = find_dopant(syms[i]);
      out << "    { \"symbol\": \"" << json_escape(syms[i]) << "\", \"name\": \""
          << json_escape(d->name) << "\", \"type\": \"" << dop_type_name(d->type)
          << "\", \"in_net_doping\": " << (d->type != DopType::neutral ? "true" : "false")
          << " }";
      out << (i + 1 < syms.size() ? ",\n" : "\n");
    }
  }
  out << "  ],\n";

  out << "  \"net_doping\": {\n";
  out << "    \"field\": \"NetDoping\",\n";
  out << "    \"convention\": \"ND-NA\",\n";
  out << "    \"positive_means\": \"n-type\",\n";
  out << "    \"uses_active_concentration\": true,\n";
  out << "    \"excluded_types\": [\"neutral\"]\n";
  out << "  },\n";

  out << "  \"fields\": {\n";
  out << "    \"cell\": [";
  for (std::size_t i = 0; i < cell_field_names.size(); ++i) {
    out << "\"" << json_escape(cell_field_names[i]) << "\"";
    if (i + 1 < cell_field_names.size()) out << ", ";
  }
  out << "],\n";
  out << "    \"point\": [";
  for (std::size_t i = 0; i < point_field_names.size(); ++i) {
    out << "\"" << json_escape(point_field_names[i]) << "\"";
    if (i + 1 < point_field_names.size()) out << ", ";
  }
  out << "]\n";
  out << "  },\n";

  std::snprintf(buf, sizeof(buf), "%.17g", st.last_temp);
  out << "  \"last_temp_k\": " << buf << "\n";
  out << "}\n";

  if (!out) throw std::runtime_error("write_device: write failed: " + meta_path);
}

}  // namespace cp
