#include <cstdio>
#include <fstream>
#include <string>

#include "cprocess/gmsh_reader.hpp"
#include "test_util.hpp"

using namespace cp;

// Single tet (unit corner tet, volume 1/6) with one tagged surface "top",
// in both MSH 2.2 and MSH 4.1 ASCII.
static const char* kMsh22 = R"(
$MeshFormat
2.2 0 8
$EndMeshFormat
$PhysicalNames
2
2 7 "top"
3 1 "sub"
$EndPhysicalNames
$Nodes
4
1 0 0 0
2 1 0 0
3 0 1 0
4 0 0 1
$EndNodes
$Elements
2
1 2 2 7 1 1 2 3
2 4 2 1 1 1 2 3 4
$EndElements
)";

static const char* kMsh41 = R"(
$MeshFormat
4.1 0 8
$EndMeshFormat
$PhysicalNames
2
2 7 "top"
3 1 "sub"
$EndPhysicalNames
$Entities
0 0 1 1
1 0 0 0 1 1 0 1 7 0
1 0 0 0 1 1 1 1 1 1 1
$EndEntities
$Nodes
1 4 1 4
2 1 0 4
1
2
3
4
0 0 0
1 0 0
0 1 0
0 0 1
$EndNodes
$Elements
2 2 1 2
2 1 2 1
1 1 2 3
3 1 4 1
2 1 2 3 4
$EndElements
)";

static std::string write_tmp(const char* name, const char* content) {
  std::string path = std::string("./") + name;
  std::ofstream f(path);
  f << content;
  return path;
}

static void check_mesh(const Mesh& m, double scale) {
  CHECK(m.cells.size() == 1);
  CHECK(m.nodes.size() == 4);
  const double s3 = scale * scale * scale;
  CHECK_NEAR(m.cell_vol[0], s3 / 6.0, 1e-9 * s3);
  CHECK(m.cell_region[0] == 1);
  CHECK(m.region_names.at(1) == "sub");

  const int top = m.find_patch("top");
  const int dflt = m.find_patch("default");
  CHECK(top >= 0);
  CHECK(dflt >= 0);
  int ntop = 0, ndflt = 0;
  double atop = 0;
  for (const auto& f : m.faces) {
    CHECK(f.neigh < 0);  // single tet: all faces are boundary
    if (f.patch == top) { ++ntop; atop += norm(f.S); }
    if (f.patch == dflt) ++ndflt;
  }
  CHECK(ntop == 1);
  CHECK(ndflt == 3);
  CHECK_NEAR(atop, 0.5 * scale * scale, 1e-9 * scale * scale);
}

int main() {
  {
    const std::string p = write_tmp("t22.msh", kMsh22);
    Mesh m = read_gmsh(p, 1.0);
    check_mesh(m, 1.0);
    std::remove(p.c_str());
  }
  {
    const std::string p = write_tmp("t41.msh", kMsh41);
    Mesh m = read_gmsh(p, 1e-4);  // mesh drawn in um
    check_mesh(m, 1e-4);
    std::remove(p.c_str());
  }
  std::printf("gmsh tests passed\n");
  return 0;
}
