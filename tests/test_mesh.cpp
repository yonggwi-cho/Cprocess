#include <cstdio>
#include <map>

#include "cprocess/mesh.hpp"
#include "test_util.hpp"

using namespace cp;

int main() {
  const double lx = 1.2, ly = 0.8, lz = 2.0;
  const int nx = 3, ny = 4, nz = 5;
  Mesh m = make_box_mesh(0, lx, 0, ly, 0, lz, nx, ny, nz);

  CHECK(m.cells.size() == static_cast<std::size_t>(nx * ny * nz * 6));
  CHECK(m.nodes.size() == static_cast<std::size_t>((nx + 1) * (ny + 1) * (nz + 1)));

  // All volumes positive, total volume exact.
  double vol = 0;
  for (double v : m.cell_vol) {
    CHECK(v > 0);
    vol += v;
  }
  CHECK_NEAR(vol, lx * ly * lz, 1e-12 * lx * ly * lz);

  // Closed surface per cell: sum of outward face area vectors is zero.
  std::vector<Vec3> ssum(m.cells.size());
  for (const auto& f : m.faces) {
    ssum[f.owner] += f.S;
    if (f.neigh >= 0) ssum[f.neigh] -= f.S;
  }
  for (const auto& s : ssum) CHECK(norm(s) < 1e-12);

  // Boundary patches: counts and areas.
  std::map<int, double> parea;
  std::map<int, int> pcount;
  for (const auto& f : m.faces) {
    if (f.neigh >= 0) {
      CHECK(f.patch == -1);
      continue;
    }
    CHECK(f.patch >= 0 && f.patch < 6);
    parea[f.patch] += norm(f.S);
    pcount[f.patch]++;
  }
  CHECK_NEAR(parea[m.find_patch("xmin")], ly * lz, 1e-12);
  CHECK_NEAR(parea[m.find_patch("xmax")], ly * lz, 1e-12);
  CHECK_NEAR(parea[m.find_patch("ymin")], lx * lz, 1e-12);
  CHECK_NEAR(parea[m.find_patch("ymax")], lx * lz, 1e-12);
  CHECK_NEAR(parea[m.find_patch("zmin")], lx * ly, 1e-12);
  CHECK_NEAR(parea[m.find_patch("zmax")], lx * ly, 1e-12);
  CHECK(pcount[m.find_patch("xmin")] == 2 * ny * nz);
  CHECK(pcount[m.find_patch("zmax")] == 2 * nx * ny);

  // Two-point scheme requires positive orthogonality on internal faces.
  const double ortho = m.min_orthogonality();
  std::printf("min orthogonality = %.4f\n", ortho);
  CHECK(ortho > 0.05);

  // finalize() must fix an inverted tet.
  Mesh t;
  t.nodes = {{0, 0, 0}, {1, 0, 0}, {0, 1, 0}, {0, 0, 1}};
  t.cells = {{0, 1, 3, 2}};  // negative volume as given
  t.cell_region = {1};
  t.finalize();
  CHECK_NEAR(t.cell_vol[0], 1.0 / 6.0, 1e-15);

  std::printf("mesh tests passed\n");
  return 0;
}
