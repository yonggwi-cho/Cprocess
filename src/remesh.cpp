#include "cprocess/remesh.hpp"

#include <algorithm>
#include <cmath>

#include "cprocess/topology.hpp"

namespace cp {

namespace {
constexpr int kEdges[6][2] = {{0, 1}, {0, 2}, {0, 3}, {1, 2}, {1, 3}, {2, 3}};

double signed_vol(const Vec3& a, const Vec3& b, const Vec3& c, const Vec3& d) {
  return dot(cross(b - a, c - a), d - a) / 6.0;
}
}  // namespace

double tet_quality(const Vec3& a, const Vec3& b, const Vec3& c, const Vec3& d) {
  const double v = signed_vol(a, b, c, d);
  const Vec3 e[6] = {b - a, c - a, d - a, c - b, d - b, d - c};
  double l2 = 0;
  for (const auto& ev : e) l2 += dot(ev, ev);
  const double lrms = std::sqrt(l2 / 6.0);
  if (lrms <= 0) return 0.0;
  // 6*sqrt(2) normalizes a regular tet to 1.
  return 6.0 * std::sqrt(2.0) * v / (lrms * lrms * lrms);
}

QualityStats mesh_quality(const Mesh& m, double sliver_thresh) {
  QualityStats s;
  const int nc = static_cast<int>(m.cells.size());
  if (nc == 0) return s;
  double sum = 0;
  for (int ci = 0; ci < nc; ++ci) {
    const auto& c = m.cells[ci];
    const double q = tet_quality(m.nodes[c[0]], m.nodes[c[1]], m.nodes[c[2]],
                                 m.nodes[c[3]]);
    sum += q;
    if (q < s.min_q) { s.min_q = q; s.worst_cell = ci; }
    if (q < sliver_thresh) ++s.n_sliver;
  }
  s.mean_q = sum / nc;
  return s;
}

int edge_split(Mesh& m, int a, int b) {
  std::vector<int> inc;
  const int nc = static_cast<int>(m.cells.size());
  for (int ci = 0; ci < nc; ++ci) {
    const auto& c = m.cells[ci];
    bool ha = false, hb = false;
    for (int v : c) { ha |= (v == a); hb |= (v == b); }
    if (ha && hb) inc.push_back(ci);
  }
  if (inc.empty()) return -1;

  const int mid = static_cast<int>(m.nodes.size());
  m.nodes.push_back(0.5 * (m.nodes[a] + m.nodes[b]));

  for (int ci : inc) {
    std::array<int, 4> c1 = m.cells[ci];  // a-side: replace b with mid
    std::array<int, 4> c2 = m.cells[ci];  // b-side: replace a with mid
    for (auto& v : c1) if (v == b) v = mid;
    for (auto& v : c2) if (v == a) v = mid;
    m.cells[ci] = c1;
    m.cells.push_back(c2);
    m.cell_region.push_back(m.cell_region[ci]);
  }
  return mid;
}

int laplacian_smooth(Mesh& m, int iters, double omega) {
  MeshTopology topo;
  topo.build(m);
  const int nn = static_cast<int>(m.nodes.size());

  // Deduplicate adjacency so each neighbour is counted once.
  std::vector<std::vector<int>> adj(nn);
  for (int i = 0; i < nn; ++i) {
    adj[i] = topo.node_adj[i];
    std::sort(adj[i].begin(), adj[i].end());
    adj[i].erase(std::unique(adj[i].begin(), adj[i].end()), adj[i].end());
  }

  auto no_inversion = [&](int node, const Vec3& pos) {
    for (int ci : topo.node_cells[node]) {
      const auto& c = m.cells[ci];
      Vec3 p[4] = {m.nodes[c[0]], m.nodes[c[1]], m.nodes[c[2]], m.nodes[c[3]]};
      for (int v = 0; v < 4; ++v) if (c[v] == node) p[v] = pos;
      if (signed_vol(p[0], p[1], p[2], p[3]) <= 0) return false;
    }
    return true;
  };

  int moved = 0;
  for (int it = 0; it < iters; ++it) {
    for (int i = 0; i < nn; ++i) {
      if (topo.node_boundary[i] || topo.node_interface[i]) continue;
      if (adj[i].empty()) continue;
      Vec3 avg{};
      for (int j : adj[i]) avg += m.nodes[j];
      avg = (1.0 / adj[i].size()) * avg;
      const Vec3 cand = (1.0 - omega) * m.nodes[i] + omega * avg;
      if (no_inversion(i, cand)) { m.nodes[i] = cand; ++moved; }
    }
  }
  return moved;
}

}  // namespace cp
