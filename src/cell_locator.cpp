#include "cprocess/cell_locator.hpp"

#include <algorithm>
#include <cmath>

namespace cp {

namespace {

double signed_vol(const Vec3& a, const Vec3& b, const Vec3& c, const Vec3& d) {
  return dot(cross(b - a, c - a), d - a) / 6.0;
}

// Cells are positively oriented (Mesh::finalize), so p is inside iff every
// vertex-replaced sub-tet keeps a non-negative volume.
bool inside_tet(const Mesh& m, int ci, const Vec3& p) {
  const auto& c = m.cells[ci];
  const Vec3 &a = m.nodes[c[0]], &b = m.nodes[c[1]], &cc = m.nodes[c[2]],
             &d = m.nodes[c[3]];
  const double tol = -1e-9 * m.cell_vol[ci];
  return signed_vol(p, b, cc, d) >= tol && signed_vol(a, p, cc, d) >= tol &&
         signed_vol(a, b, p, d) >= tol && signed_vol(a, b, cc, p) >= tol;
}

}  // namespace

CellLocator::CellLocator(const Mesh& mesh) : mesh_(mesh) {
  const BBox bb = mesh.bbox();
  lo_ = bb.lo;
  hi_ = bb.hi;
  const double ex = std::max(hi_.x - lo_.x, 1e-300);
  const double ey = std::max(hi_.y - lo_.y, 1e-300);
  const double ez = std::max(hi_.z - lo_.z, 1e-300);

  // Roughly one cell per voxel, axes proportional to the extents.
  const double s = std::cbrt(static_cast<double>(mesh.cells.size()) /
                             (ex * ey * ez));
  auto naxis = [&](double e) {
    return std::max(1, std::min(512, static_cast<int>(std::ceil(e * s))));
  };
  nx_ = naxis(ex);
  ny_ = naxis(ey);
  nz_ = naxis(ez);
  ix_ = nx_ / ex;
  iy_ = ny_ / ey;
  iz_ = nz_ / ez;

  auto clampi = [](int v, int n) { return std::max(0, std::min(n - 1, v)); };
  const int nvox = nx_ * ny_ * nz_;
  std::vector<int> count(nvox + 1, 0);

  auto cell_range = [&](int ci, int r[6]) {
    const auto& c = mesh.cells[ci];
    Vec3 clo = mesh.nodes[c[0]], chi = clo;
    for (int k = 1; k < 4; ++k) {
      const Vec3& p = mesh.nodes[c[k]];
      clo.x = std::min(clo.x, p.x); chi.x = std::max(chi.x, p.x);
      clo.y = std::min(clo.y, p.y); chi.y = std::max(chi.y, p.y);
      clo.z = std::min(clo.z, p.z); chi.z = std::max(chi.z, p.z);
    }
    r[0] = clampi(static_cast<int>((clo.x - lo_.x) * ix_), nx_);
    r[1] = clampi(static_cast<int>((chi.x - lo_.x) * ix_), nx_);
    r[2] = clampi(static_cast<int>((clo.y - lo_.y) * iy_), ny_);
    r[3] = clampi(static_cast<int>((chi.y - lo_.y) * iy_), ny_);
    r[4] = clampi(static_cast<int>((clo.z - lo_.z) * iz_), nz_);
    r[5] = clampi(static_cast<int>((chi.z - lo_.z) * iz_), nz_);
  };

  for (std::size_t ci = 0; ci < mesh.cells.size(); ++ci) {
    int r[6];
    cell_range(static_cast<int>(ci), r);
    for (int k = r[4]; k <= r[5]; ++k)
      for (int j = r[2]; j <= r[3]; ++j)
        for (int i = r[0]; i <= r[1]; ++i)
          ++count[(k * ny_ + j) * nx_ + i + 1];
  }
  for (int v = 0; v < nvox; ++v) count[v + 1] += count[v];
  start_ = count;
  items_.resize(start_[nvox]);
  std::vector<int> fill(start_.begin(), start_.end() - 1);
  for (std::size_t ci = 0; ci < mesh.cells.size(); ++ci) {
    int r[6];
    cell_range(static_cast<int>(ci), r);
    for (int k = r[4]; k <= r[5]; ++k)
      for (int j = r[2]; j <= r[3]; ++j)
        for (int i = r[0]; i <= r[1]; ++i)
          items_[fill[(k * ny_ + j) * nx_ + i]++] = static_cast<int>(ci);
  }
}

int CellLocator::locate(const Vec3& p) const {
  if (p.x < lo_.x || p.x > hi_.x || p.y < lo_.y || p.y > hi_.y ||
      p.z < lo_.z || p.z > hi_.z)
    return -1;
  auto clampi = [](int v, int n) { return std::max(0, std::min(n - 1, v)); };
  const int i = clampi(static_cast<int>((p.x - lo_.x) * ix_), nx_);
  const int j = clampi(static_cast<int>((p.y - lo_.y) * iy_), ny_);
  const int k = clampi(static_cast<int>((p.z - lo_.z) * iz_), nz_);
  const int v = (k * ny_ + j) * nx_ + i;
  for (int s = start_[v]; s < start_[v + 1]; ++s)
    if (inside_tet(mesh_, items_[s], p)) return items_[s];
  return -1;
}

}  // namespace cp
