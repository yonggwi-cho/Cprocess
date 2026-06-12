#include "cprocess/mesh.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace cp {

static double signed_volume(const Vec3& a, const Vec3& b, const Vec3& c,
                            const Vec3& d) {
  return dot(cross(b - a, c - a), d - a) / 6.0;
}

void Mesh::finalize() {
  const int nc = static_cast<int>(cells.size());
  if (nodes.empty() || nc == 0) throw std::runtime_error("mesh: empty mesh");
  if (static_cast<int>(cell_region.size()) != nc) cell_region.assign(nc, 0);

  // Orient every tet to positive volume.
  for (auto& c : cells) {
    const double v = signed_volume(nodes[c[0]], nodes[c[1]], nodes[c[2]], nodes[c[3]]);
    if (v < 0) std::swap(c[2], c[3]);
  }

  // Build faces. For a positively oriented tet (a,b,c,d) the outward faces
  // are (a,c,b), (a,b,d), (a,d,c), (b,c,d).
  static const int kFV[4][3] = {{0, 2, 1}, {0, 1, 3}, {0, 3, 2}, {1, 2, 3}};
  faces.clear();
  face_lookup.clear();
  face_lookup.reserve(cells.size() * 2 + 16);
  for (int ci = 0; ci < nc; ++ci) {
    for (const auto& fv : kFV) {
      std::array<int, 3> tri = {cells[ci][fv[0]], cells[ci][fv[1]], cells[ci][fv[2]]};
      std::array<int, 3> key = tri;
      std::sort(key.begin(), key.end());
      auto it = face_lookup.find(key);
      if (it == face_lookup.end()) {
        Face f;
        f.n = tri;
        f.owner = ci;
        face_lookup.emplace(key, static_cast<int>(faces.size()));
        faces.push_back(f);
      } else {
        Face& f = faces[it->second];
        if (f.neigh >= 0)
          throw std::runtime_error("mesh: non-manifold face (shared by >2 cells)");
        f.neigh = ci;
      }
    }
  }

  // Geometry.
  cell_vol.resize(nc);
  cell_cent.resize(nc);
  for (int ci = 0; ci < nc; ++ci) {
    const auto& c = cells[ci];
    const Vec3 &a = nodes[c[0]], &b = nodes[c[1]], &p = nodes[c[2]], &q = nodes[c[3]];
    cell_vol[ci] = signed_volume(a, b, p, q);
    if (cell_vol[ci] <= 0.0)
      throw std::runtime_error("mesh: degenerate cell (zero volume)");
    cell_cent[ci] = 0.25 * (a + b + p + q);
  }
  for (auto& f : faces) {
    const Vec3 &a = nodes[f.n[0]], &b = nodes[f.n[1]], &c = nodes[f.n[2]];
    f.S = 0.5 * cross(b - a, c - a);
    f.c = (1.0 / 3.0) * (a + b + c);
  }
}

int Mesh::add_patch(const std::string& name) {
  const int id = find_patch(name);
  if (id >= 0) return id;
  patch_names.push_back(name);
  return static_cast<int>(patch_names.size()) - 1;
}

int Mesh::find_patch(const std::string& name) const {
  for (std::size_t i = 0; i < patch_names.size(); ++i)
    if (patch_names[i] == name) return static_cast<int>(i);
  return -1;
}

BBox Mesh::bbox() const {
  BBox b;
  const double inf = std::numeric_limits<double>::infinity();
  b.lo = {inf, inf, inf};
  b.hi = {-inf, -inf, -inf};
  for (const auto& p : nodes) {
    b.lo.x = std::min(b.lo.x, p.x); b.hi.x = std::max(b.hi.x, p.x);
    b.lo.y = std::min(b.lo.y, p.y); b.hi.y = std::max(b.hi.y, p.y);
    b.lo.z = std::min(b.lo.z, p.z); b.hi.z = std::max(b.hi.z, p.z);
  }
  return b;
}

double Mesh::total_volume() const {
  double v = 0;
  for (double cv : cell_vol) v += cv;
  return v;
}

double Mesh::min_orthogonality() const {
  double m = 1.0;
  for (const auto& f : faces) {
    if (f.neigh < 0) continue;
    const Vec3 d = cell_cent[f.neigh] - cell_cent[f.owner];
    const double den = norm(f.S) * norm(d);
    if (den > 0) m = std::min(m, dot(f.S, d) / den);
  }
  return m;
}

std::vector<int> Mesh::region_tags() const {
  std::vector<int> t(cell_region);
  std::sort(t.begin(), t.end());
  t.erase(std::unique(t.begin(), t.end()), t.end());
  return t;
}

}  // namespace cp
