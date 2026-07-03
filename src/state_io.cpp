#include "cprocess/state_io.hpp"

#include <cstdint>
#include <fstream>
#include <stdexcept>

namespace cp {

namespace {

void write_i64(std::ofstream& f, std::int64_t v) {
  f.write(reinterpret_cast<const char*>(&v), sizeof(v));
}
void write_f64(std::ofstream& f, double v) {
  f.write(reinterpret_cast<const char*>(&v), sizeof(v));
}
void write_str(std::ofstream& f, const std::string& s) {
  write_i64(f, static_cast<std::int64_t>(s.size()));
  if (!s.empty()) f.write(s.data(), static_cast<std::streamsize>(s.size()));
}

std::int64_t read_i64(std::ifstream& f) {
  std::int64_t v = 0;
  f.read(reinterpret_cast<char*>(&v), sizeof(v));
  if (!f) throw std::runtime_error("not a CPRC1 state file");
  return v;
}
double read_f64(std::ifstream& f) {
  double v = 0;
  f.read(reinterpret_cast<char*>(&v), sizeof(v));
  if (!f) throw std::runtime_error("not a CPRC1 state file");
  return v;
}
std::string read_str(std::ifstream& f) {
  std::int64_t n = read_i64(f);
  std::string s;
  if (n > 0) {
    s.resize(static_cast<std::size_t>(n));
    f.read(&s[0], n);
    if (!f) throw std::runtime_error("not a CPRC1 state file");
  }
  return s;
}

}  // namespace

void write_state(const SimState& st, const std::string& path) {
  std::ofstream f(path, std::ios::binary);
  if (!f) throw std::runtime_error("cannot write '" + path + "'");

  f.write("CPRC1", 5);
  write_i64(f, 1);

  const auto& nodes = st.mesh.nodes;
  write_i64(f, static_cast<std::int64_t>(nodes.size()));
  for (const auto& n : nodes) {
    write_f64(f, n.x);
    write_f64(f, n.y);
    write_f64(f, n.z);
  }

  const auto& cells = st.mesh.cells;
  const std::size_t nc = cells.size();
  write_i64(f, static_cast<std::int64_t>(nc));
  for (const auto& c : cells) {
    for (int id : c) write_i64(f, id);
  }

  for (int r : st.mesh.cell_region) write_i64(f, r);

  write_i64(f, static_cast<std::int64_t>(st.mesh.region_names.size()));
  for (const auto& [tag, name] : st.mesh.region_names) {
    write_i64(f, tag);
    write_str(f, name);
  }

  write_i64(f, static_cast<std::int64_t>(st.mesh.patch_names.size()));
  for (const auto& name : st.mesh.patch_names) write_str(f, name);

  write_i64(f, static_cast<std::int64_t>(st.region_material.size()));
  for (const auto& [tag, mat] : st.region_material) {
    write_i64(f, tag);
    write_str(f, mat);
  }

  write_i64(f, static_cast<std::int64_t>(st.fields.size()));
  for (const auto& [sym, conc] : st.fields) {
    write_str(f, sym);
    if (conc.size() == nc) {
      for (double v : conc) write_f64(f, v);
    } else {
      for (std::size_t i = 0; i < nc; ++i) write_f64(f, 0.0);
    }
  }

  write_i64(f, static_cast<std::int64_t>(st.bcs.size()));
  for (const auto& bc : st.bcs) {
    write_str(f, bc.species);
    write_i64(f, bc.patch);
    write_f64(f, bc.conc);
  }

  write_f64(f, st.last_temp);

  if (!f) throw std::runtime_error("cannot write '" + path + "'");
}

void read_state(SimState& st, const std::string& path) {
  std::ifstream f(path, std::ios::binary);
  if (!f) throw std::runtime_error("cannot read '" + path + "'");

  char magic[5] = {0};
  f.read(magic, 5);
  if (!f || std::string(magic, 5) != "CPRC1")
    throw std::runtime_error("not a CPRC1 state file");

  std::int64_t version = read_i64(f);
  if (version != 1) throw std::runtime_error("unsupported state version");

  SimState out;

  std::int64_t n_nodes = read_i64(f);
  out.mesh.nodes.resize(static_cast<std::size_t>(n_nodes));
  for (auto& n : out.mesh.nodes) {
    n.x = read_f64(f);
    n.y = read_f64(f);
    n.z = read_f64(f);
  }

  std::int64_t n_cells = read_i64(f);
  out.mesh.cells.resize(static_cast<std::size_t>(n_cells));
  for (auto& c : out.mesh.cells) {
    for (int& id : c) id = static_cast<int>(read_i64(f));
  }

  out.mesh.cell_region.resize(static_cast<std::size_t>(n_cells));
  for (auto& r : out.mesh.cell_region) r = static_cast<int>(read_i64(f));

  std::int64_t n_region_names = read_i64(f);
  for (std::int64_t i = 0; i < n_region_names; ++i) {
    int tag = static_cast<int>(read_i64(f));
    std::string name = read_str(f);
    out.mesh.region_names[tag] = name;
  }

  std::int64_t n_patch_names = read_i64(f);
  out.mesh.patch_names.resize(static_cast<std::size_t>(n_patch_names));
  for (auto& name : out.mesh.patch_names) name = read_str(f);

  std::int64_t n_region_material = read_i64(f);
  for (std::int64_t i = 0; i < n_region_material; ++i) {
    int tag = static_cast<int>(read_i64(f));
    std::string mat = read_str(f);
    out.region_material[tag] = mat;
  }

  std::int64_t n_fields = read_i64(f);
  for (std::int64_t i = 0; i < n_fields; ++i) {
    std::string sym = read_str(f);
    std::vector<double> conc(static_cast<std::size_t>(n_cells));
    for (auto& v : conc) v = read_f64(f);
    out.fields[sym] = std::move(conc);
  }

  std::int64_t n_bcs = read_i64(f);
  out.bcs.resize(static_cast<std::size_t>(n_bcs));
  for (auto& bc : out.bcs) {
    bc.species = read_str(f);
    bc.patch = static_cast<int>(read_i64(f));
    bc.conc = read_f64(f);
  }

  out.last_temp = read_f64(f);

  if (!f) throw std::runtime_error("not a CPRC1 state file");

  out.mesh.finalize();
  out.has_mesh = true;
  out.has_stack = false;
  out.stack = Mesh();
  out.stack_cell_mat.clear();
  out.mat_table.clear();
  out.stack_resist_z0 = 0;
  out.layer_stack.clear();

  st = std::move(out);
}

}  // namespace cp
