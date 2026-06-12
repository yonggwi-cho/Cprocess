#include "cprocess/gmsh_reader.hpp"

#include <algorithm>
#include <fstream>
#include <map>
#include <ostream>
#include <sstream>
#include <stdexcept>

namespace cp {

namespace {

struct Reader {
  std::ifstream in;
  int lineno = 0;
  explicit Reader(const std::string& path) : in(path) {}

  // Next non-empty line (trailing \r stripped). Returns false on EOF.
  bool next(std::string& line) {
    while (std::getline(in, line)) {
      ++lineno;
      while (!line.empty() && (line.back() == '\r' || line.back() == '\n'))
        line.pop_back();
      if (!line.empty()) return true;
    }
    return false;
  }
  [[noreturn]] void fail(const std::string& msg) const {
    throw std::runtime_error("gmsh: line " + std::to_string(lineno) + ": " + msg);
  }
};

struct Tri { int phys; int n[3]; };
struct Tet { int phys; int n[4]; };

int first_order_nodes(int etype) {
  switch (etype) {
    case 1: return 2;   // line
    case 2: return 3;   // triangle
    case 3: return 4;   // quad
    case 4: return 4;   // tetrahedron
    case 15: return 1;  // point
    default: return -1;
  }
}

}  // namespace

Mesh read_gmsh(const std::string& path, double scale, std::ostream* log) {
  Reader r(path);
  if (!r.in) throw std::runtime_error("gmsh: cannot open file: " + path);

  double version = 0;
  std::map<int, int> node_index;                 // gmsh node tag -> index
  std::vector<Vec3> nodes;
  std::vector<Tri> tris;
  std::vector<Tet> tets;
  std::map<int, std::string> phys_names[4];      // by dimension
  std::map<std::pair<int, int>, int> ent_phys;   // (dim, entity tag) -> phys tag
  int skipped_higher = 0;

  std::string line;
  while (r.next(line)) {
    if (line[0] != '$') continue;
    const std::string section = line.substr(1);

    auto skip_section = [&]() {
      while (r.next(line))
        if (line == "$End" + section) return;
      r.fail("missing $End" + section);
    };

    if (section == "MeshFormat") {
      if (!r.next(line)) r.fail("truncated $MeshFormat");
      std::istringstream is(line);
      int ftype = -1, dsize = 0;
      is >> version >> ftype >> dsize;
      if (ftype != 0) r.fail("binary .msh not supported; export as ASCII");
      if (!((version >= 2.0 && version < 3.0) || (version >= 4.0 && version < 5.0)))
        r.fail("unsupported MSH version (use 2.2 or 4.1 ASCII)");
      skip_section();
    } else if (section == "PhysicalNames") {
      if (!r.next(line)) r.fail("truncated $PhysicalNames");
      const int n = std::stoi(line);
      for (int i = 0; i < n; ++i) {
        if (!r.next(line)) r.fail("truncated $PhysicalNames");
        std::istringstream is(line);
        int dim = 0, tag = 0;
        is >> dim >> tag;
        const auto q1 = line.find('"');
        const auto q2 = line.rfind('"');
        std::string name = (q1 != std::string::npos && q2 > q1)
                               ? line.substr(q1 + 1, q2 - q1 - 1)
                               : std::string();
        if (dim >= 0 && dim <= 3) phys_names[dim][tag] = name;
      }
      skip_section();
    } else if (section == "Entities" && version >= 4.0) {
      if (!r.next(line)) r.fail("truncated $Entities");
      long np = 0, ncv = 0, ns = 0, nv = 0;
      { std::istringstream is(line); is >> np >> ncv >> ns >> nv; }
      for (long i = 0; i < np + ncv; ++i)
        if (!r.next(line)) r.fail("truncated $Entities");
      auto read_ent = [&](int dim) {
        std::istringstream is(line);
        int tag; double bb;
        is >> tag;
        for (int k = 0; k < 6; ++k) is >> bb;
        int nphys = 0;
        is >> nphys;
        int phys = 0;
        for (int k = 0; k < nphys; ++k) {
          int p; is >> p;
          if (k == 0) phys = p;
        }
        ent_phys[{dim, tag}] = phys;
      };
      for (long i = 0; i < ns; ++i) {
        if (!r.next(line)) r.fail("truncated $Entities");
        read_ent(2);
      }
      for (long i = 0; i < nv; ++i) {
        if (!r.next(line)) r.fail("truncated $Entities");
        read_ent(3);
      }
      skip_section();
    } else if (section == "Nodes") {
      if (!r.next(line)) r.fail("truncated $Nodes");
      if (version < 3.0) {
        const long n = std::stol(line);
        for (long i = 0; i < n; ++i) {
          if (!r.next(line)) r.fail("truncated $Nodes");
          std::istringstream is(line);
          long tag; double x, y, z;
          if (!(is >> tag >> x >> y >> z)) r.fail("bad node line");
          node_index[static_cast<int>(tag)] = static_cast<int>(nodes.size());
          nodes.push_back({x * scale, y * scale, z * scale});
        }
      } else {
        long nblocks, nnodes, mn, mx;
        { std::istringstream is(line); is >> nblocks >> nnodes >> mn >> mx; }
        for (long b = 0; b < nblocks; ++b) {
          if (!r.next(line)) r.fail("truncated $Nodes");
          int dim, etag, param;
          long nb;
          { std::istringstream is(line); is >> dim >> etag >> param >> nb; }
          if (param != 0) r.fail("parametric nodes not supported");
          std::vector<long> tags(nb);
          for (long i = 0; i < nb; ++i) {
            if (!r.next(line)) r.fail("truncated $Nodes");
            tags[i] = std::stol(line);
          }
          for (long i = 0; i < nb; ++i) {
            if (!r.next(line)) r.fail("truncated $Nodes");
            std::istringstream is(line);
            double x, y, z;
            if (!(is >> x >> y >> z)) r.fail("bad node coordinates");
            node_index[static_cast<int>(tags[i])] = static_cast<int>(nodes.size());
            nodes.push_back({x * scale, y * scale, z * scale});
          }
        }
      }
      skip_section();
    } else if (section == "Elements") {
      if (!r.next(line)) r.fail("truncated $Elements");
      if (version < 3.0) {
        const long n = std::stol(line);
        for (long i = 0; i < n; ++i) {
          if (!r.next(line)) r.fail("truncated $Elements");
          std::istringstream is(line);
          long id; int etype, ntags;
          if (!(is >> id >> etype >> ntags)) r.fail("bad element line");
          int phys = 0, tagv;
          for (int k = 0; k < ntags; ++k) {
            is >> tagv;
            if (k == 0) phys = tagv;
          }
          if (etype == 2) {
            Tri t{phys, {0, 0, 0}};
            if (!(is >> t.n[0] >> t.n[1] >> t.n[2])) r.fail("bad triangle");
            tris.push_back(t);
          } else if (etype == 4) {
            Tet t{phys, {0, 0, 0, 0}};
            if (!(is >> t.n[0] >> t.n[1] >> t.n[2] >> t.n[3])) r.fail("bad tet");
            tets.push_back(t);
          } else if (first_order_nodes(etype) < 0) {
            ++skipped_higher;
          }
        }
      } else {
        long nblocks, nelems, mn, mx;
        { std::istringstream is(line); is >> nblocks >> nelems >> mn >> mx; }
        for (long b = 0; b < nblocks; ++b) {
          if (!r.next(line)) r.fail("truncated $Elements");
          int dim, etag, etype;
          long ne;
          { std::istringstream is(line); is >> dim >> etag >> etype >> ne; }
          int phys = 0;
          auto it = ent_phys.find({dim, etag});
          if (it != ent_phys.end()) phys = it->second;
          for (long i = 0; i < ne; ++i) {
            if (!r.next(line)) r.fail("truncated $Elements");
            std::istringstream is(line);
            long id;
            is >> id;
            if (etype == 2) {
              Tri t{phys, {0, 0, 0}};
              if (!(is >> t.n[0] >> t.n[1] >> t.n[2])) r.fail("bad triangle");
              tris.push_back(t);
            } else if (etype == 4) {
              Tet t{phys != 0 ? phys : etag, {0, 0, 0, 0}};
              if (!(is >> t.n[0] >> t.n[1] >> t.n[2] >> t.n[3])) r.fail("bad tet");
              tets.push_back(t);
            } else if (first_order_nodes(etype) < 0) {
              ++skipped_higher;
            }
          }
        }
      }
      skip_section();
    } else {
      skip_section();
    }
  }

  if (skipped_higher > 0)
    throw std::runtime_error(
        "gmsh: higher-order elements found; generate a 1st-order mesh "
        "(gmsh -order 1)");
  if (tets.empty()) throw std::runtime_error("gmsh: no tetrahedra in file");

  Mesh m;
  m.nodes = std::move(nodes);
  m.cells.reserve(tets.size());
  m.cell_region.reserve(tets.size());
  auto idx = [&](int tag) {
    auto it = node_index.find(tag);
    if (it == node_index.end())
      throw std::runtime_error("gmsh: element references unknown node " +
                               std::to_string(tag));
    return it->second;
  };
  for (const auto& t : tets) {
    m.cells.push_back({idx(t.n[0]), idx(t.n[1]), idx(t.n[2]), idx(t.n[3])});
    m.cell_region.push_back(t.phys);
  }
  m.region_names = phys_names[3];
  m.finalize();

  // Assign boundary patches from tagged surface triangles.
  std::map<int, int> patch_of_phys;
  int internal_tagged = 0, unmatched = 0;
  for (const auto& t : tris) {
    if (t.phys == 0) continue;
    std::array<int, 3> key = {idx(t.n[0]), idx(t.n[1]), idx(t.n[2])};
    std::sort(key.begin(), key.end());
    auto it = m.face_lookup.find(key);
    if (it == m.face_lookup.end()) {
      ++unmatched;
      continue;
    }
    Face& f = m.faces[it->second];
    if (f.neigh >= 0) {
      ++internal_tagged;  // interface inside the volume: not a BC patch
      continue;
    }
    auto pit = patch_of_phys.find(t.phys);
    if (pit == patch_of_phys.end()) {
      std::string name;
      auto nit = phys_names[2].find(t.phys);
      name = (nit != phys_names[2].end() && !nit->second.empty())
                 ? nit->second
                 : "surface_" + std::to_string(t.phys);
      pit = patch_of_phys.emplace(t.phys, m.add_patch(name)).first;
    }
    f.patch = pit->second;
  }
  bool any_default = false;
  for (auto& f : m.faces)
    if (f.neigh < 0 && f.patch < 0) any_default = true;
  if (any_default) {
    const int dp = m.add_patch("default");
    for (auto& f : m.faces)
      if (f.neigh < 0 && f.patch < 0) f.patch = dp;
  }

  if (log) {
    *log << "[mesh] gmsh '" << path << "' (v" << version << "): "
         << m.nodes.size() << " nodes, " << m.cells.size() << " tets, "
         << m.patch_names.size() << " patches";
    if (internal_tagged) *log << ", " << internal_tagged << " internal tagged tris ignored";
    if (unmatched) *log << ", " << unmatched << " unmatched tris ignored";
    *log << "\n";
  }
  return m;
}

}  // namespace cp
