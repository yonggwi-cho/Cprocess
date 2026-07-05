#include "cprocess/ale_mover.hpp"

#include <algorithm>
#include <cmath>

namespace cp {

namespace {
double signed_vol(const Vec3& a, const Vec3& b, const Vec3& c, const Vec3& d) {
  return dot(cross(b - a, c - a), d - a) / 6.0;
}
}  // namespace

std::vector<Vec3> compute_node_normals(const Mesh& m, const MeshTopology& topo) {
  const int nn = static_cast<int>(m.nodes.size());
  std::vector<Vec3> normals(nn, {0, 0, 0});

  for (const auto& f : m.faces) {
    if (f.neigh >= 0) continue;  // internal face
    // Outward area vector S already points outward for boundary faces.
    const double area = norm(f.S);
    if (area <= 0) continue;
    const Vec3 n = (1.0 / area) * f.S;
    for (int ni : f.n) {
      if (topo.node_boundary[ni] || topo.node_interface[ni]) {
        normals[ni] += n * area;
      }
    }
  }

  // Also account for interface faces (between regions).
  for (const auto& f : m.faces) {
    if (f.neigh < 0) continue;
    const int ra = m.cell_region[f.owner], rb = m.cell_region[f.neigh];
    if (ra == rb) continue;
    const double area = norm(f.S);
    if (area <= 0) continue;
    const Vec3 n = (1.0 / area) * f.S;
    for (int ni : f.n) {
      if (topo.node_interface[ni]) normals[ni] += n * area;
    }
  }

  // Normalize.
  for (int i = 0; i < nn; ++i) {
    const double len = norm(normals[i]);
    if (len > 0) normals[i] = (1.0 / len) * normals[i];
  }
  return normals;
}

namespace {

// Shared implementation: `disp_of(i)` gives the displacement to apply to
// boundary/interface node i (the scalar and per-node overloads differ only
// in how this displacement is computed).
AleResult ale_move_impl(Mesh& m, const MeshTopology& topo,
                        const std::function<Vec3(int)>& disp_of,
                        std::function<bool(int)> mask_fn, int smooth_iters) {
  const int nn = static_cast<int>(m.nodes.size());

  auto no_inversion = [&](int node, const Vec3& pos) {
    for (int ci : topo.node_cells[node]) {
      const auto& c = m.cells[ci];
      Vec3 p[4] = {m.nodes[c[0]], m.nodes[c[1]], m.nodes[c[2]], m.nodes[c[3]]};
      for (int v = 0; v < 4; ++v) if (c[v] == node) p[v] = pos;
      if (signed_vol(p[0], p[1], p[2], p[3]) <= 0) return false;
    }
    return true;
  };

  AleResult res;

  // Move boundary/interface nodes.
  for (int i = 0; i < nn; ++i) {
    if (!topo.node_boundary[i] && !topo.node_interface[i]) continue;
    if (mask_fn && mask_fn(i)) continue;

    const Vec3 disp = disp_of(i);
    if (norm(disp) < 1e-300) continue;
    const Vec3 cand = m.nodes[i] + disp;
    if (no_inversion(i, cand)) {
      m.nodes[i] = cand;
      ++res.n_moved;
    } else {
      ++res.n_skipped;
    }
  }

  // Laplacian smooth interior nodes to propagate displacement inward.
  for (int it = 0; it < smooth_iters; ++it) {
    for (int i = 0; i < nn; ++i) {
      if (topo.node_boundary[i] || topo.node_interface[i]) continue;
      const auto& adj = topo.node_adj[i];
      if (adj.empty()) continue;
      Vec3 avg{};
      for (int j : adj) avg += m.nodes[j];
      avg = (1.0 / adj.size()) * avg;
      if (no_inversion(i, avg)) m.nodes[i] = avg;
    }
  }

  return res;
}

}  // namespace

AleResult ale_move(Mesh& m, const MeshTopology& topo,
                   const std::vector<Vec3>& node_normals,
                   double v_n, double dt,
                   std::function<bool(int)> mask_fn,
                   int smooth_iters) {
  return ale_move_impl(
      m, topo,
      [&](int i) { return (norm(node_normals[i]) < 1e-15) ? Vec3{} : (v_n * dt) * node_normals[i]; },
      std::move(mask_fn), smooth_iters);
}

AleResult ale_move(Mesh& m, const MeshTopology& topo,
                   const std::vector<Vec3>& node_normals,
                   const std::vector<double>& vn_node, double dt,
                   std::function<bool(int)> mask_fn,
                   int smooth_iters) {
  return ale_move_impl(
      m, topo,
      [&](int i) {
        return (norm(node_normals[i]) < 1e-15) ? Vec3{} : (vn_node[i] * dt) * node_normals[i];
      },
      std::move(mask_fn), smooth_iters);
}

void rescale_fields_for_volume_change(
    std::vector<std::vector<double>>& fields,
    const std::vector<double>& old_vol,
    const std::vector<double>& new_vol) {
  const int nc = static_cast<int>(old_vol.size());
  for (auto& field : fields) {
    for (int i = 0; i < nc; ++i) {
      if (new_vol[i] > 0 && old_vol[i] > 0) {
        field[i] *= old_vol[i] / new_vol[i];
      }
    }
  }
}

}  // namespace cp
