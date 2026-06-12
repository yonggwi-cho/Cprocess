#include <algorithm>
#include <cmath>
#include <stdexcept>

#include "cprocess/mesh.hpp"

namespace cp {

// Conforming Kuhn subdivision: each hex is cut into the 6 tets
// {v0, v0+e_p1, v0+e_p1+e_p2, v7} over the permutations of (x,y,z); all
// tets share the main diagonal v0-v7, so adjacent hexes match on the
// shared quad diagonals.
Mesh make_box_mesh(double x0, double x1, double y0, double y1, double z0,
                   double z1, int nx, int ny, int nz) {
  if (nx < 1 || ny < 1 || nz < 1)
    throw std::runtime_error("box mesh: nx, ny, nz must be >= 1");
  if (!(x1 > x0) || !(y1 > y0) || !(z1 > z0))
    throw std::runtime_error("box mesh: max coordinates must exceed min");

  Mesh m;
  const int px = nx + 1, py = ny + 1, pz = nz + 1;
  m.nodes.reserve(static_cast<std::size_t>(px) * py * pz);
  for (int k = 0; k < pz; ++k)
    for (int j = 0; j < py; ++j)
      for (int i = 0; i < px; ++i)
        m.nodes.push_back({x0 + (x1 - x0) * i / nx, y0 + (y1 - y0) * j / ny,
                           z0 + (z1 - z0) * k / nz});
  auto nid = [&](int i, int j, int k) { return (k * py + j) * px + i; };

  // Hex corner bit order: b = di + 2*dj + 4*dk.
  static const int kTets[6][4] = {{0, 1, 3, 7}, {0, 1, 5, 7}, {0, 2, 3, 7},
                                  {0, 2, 6, 7}, {0, 4, 5, 7}, {0, 4, 6, 7}};
  m.cells.reserve(static_cast<std::size_t>(nx) * ny * nz * 6);
  for (int k = 0; k < nz; ++k)
    for (int j = 0; j < ny; ++j)
      for (int i = 0; i < nx; ++i) {
        int v[8];
        for (int b = 0; b < 8; ++b)
          v[b] = nid(i + (b & 1), j + ((b >> 1) & 1), k + ((b >> 2) & 1));
        for (const auto& t : kTets)
          m.cells.push_back({v[t[0]], v[t[1]], v[t[2]], v[t[3]]});
      }
  m.cell_region.assign(m.cells.size(), 1);
  m.region_names[1] = "sub";
  m.finalize();

  // Classify boundary faces onto the six box planes.
  m.patch_names = {"xmin", "xmax", "ymin", "ymax", "zmin", "zmax"};
  const double span = std::max({x1 - x0, y1 - y0, z1 - z0});
  const double tol = 1e-9 * span;
  for (auto& f : m.faces) {
    if (f.neigh >= 0) continue;
    bool on[6] = {true, true, true, true, true, true};
    for (int nn : f.n) {
      const Vec3& p = m.nodes[nn];
      on[0] &= std::fabs(p.x - x0) < tol; on[1] &= std::fabs(p.x - x1) < tol;
      on[2] &= std::fabs(p.y - y0) < tol; on[3] &= std::fabs(p.y - y1) < tol;
      on[4] &= std::fabs(p.z - z0) < tol; on[5] &= std::fabs(p.z - z1) < tol;
    }
    f.patch = -1;
    for (int s = 0; s < 6; ++s)
      if (on[s]) { f.patch = s; break; }
    if (f.patch < 0)
      throw std::runtime_error("box mesh: boundary face not on any box plane");
  }
  return m;
}

}  // namespace cp
