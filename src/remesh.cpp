#include "cprocess/remesh.hpp"

#include <algorithm>
#include <cmath>
#include <numeric>
#include <set>

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

SplitResult split_edges(Mesh& m, const std::vector<std::pair<int, int>>& edges) {
  SplitResult res;
  const int nc0 = static_cast<int>(m.cells.size());
  res.cell_parent.resize(nc0);
  std::iota(res.cell_parent.begin(), res.cell_parent.end(), 0);

  MeshTopology topo;
  topo.build(m);

  std::vector<char> cell_touched(nc0, 0);

  for (const auto& [a, b] : edges) {
    const std::vector<int>& inc = topo.edge_incident(a, b);
    if (inc.empty()) { ++res.n_skipped; continue; }

    bool conflict = false;
    for (int ci : inc) if (cell_touched[ci]) { conflict = true; break; }
    if (conflict) { ++res.n_skipped; continue; }

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
      res.cell_parent.push_back(res.cell_parent[ci]);

      cell_touched[ci] = 1;
      // ci is still < nc0 sized cell_touched; new cell has no entry, but its
      // parent ci is already marked touched, which is what matters for
      // conflict checks against `inc` (all indices < nc0 here since incident
      // cells come from the topology built before any split in this pass).
    }
    ++res.n_split;
    res.last_new_node = mid;
  }

  m.finalize();
  return res;
}

int edge_split(Mesh& m, int a, int b) {
  const std::size_t n0 = m.nodes.size();
  SplitResult res = split_edges(m, {{a, b}});
  if (m.nodes.size() == n0) return -1;
  return res.last_new_node;
}

std::vector<double> redistribute_field(const std::vector<double>& conc,
                                       const std::vector<int>& cell_parent) {
  std::vector<double> out(cell_parent.size());
  for (std::size_t i = 0; i < cell_parent.size(); ++i)
    out[i] = conc[cell_parent[i]];
  return out;
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

RepairResult repair_quality(Mesh& m, double q_thresh, int max_rounds,
                            int smooth_iters) {
  RepairResult res;
  res.min_q_before = mesh_quality(m, q_thresh).min_q;

  const int nc0 = static_cast<int>(m.cells.size());
  res.cell_parent.resize(nc0);
  std::iota(res.cell_parent.begin(), res.cell_parent.end(), 0);

  QualityStats qs = mesh_quality(m, q_thresh);

  for (int round = 0; round < max_rounds; ++round) {
    if (qs.n_sliver == 0) break;

    res.n_smoothed += laplacian_smooth(m, smooth_iters, 0.5);
    m.finalize();

    qs = mesh_quality(m, q_thresh);
    if (qs.n_sliver == 0) break;

    // Collect the longest edge of every sliver cell, deduplicated.
    std::vector<std::pair<int, int>> candidates;
    std::set<std::pair<int, int>> seen;
    const int nc = static_cast<int>(m.cells.size());
    for (int ci = 0; ci < nc; ++ci) {
      const auto& c = m.cells[ci];
      const double q = tet_quality(m.nodes[c[0]], m.nodes[c[1]], m.nodes[c[2]],
                                   m.nodes[c[3]]);
      if (q >= q_thresh) continue;
      int best = -1;
      double best_len2 = -1;
      for (int e = 0; e < 6; ++e) {
        const int va = c[kEdges[e][0]], vb = c[kEdges[e][1]];
        const Vec3 d = m.nodes[va] - m.nodes[vb];
        const double len2 = dot(d, d);
        if (len2 > best_len2) { best_len2 = len2; best = e; }
      }
      int a = c[kEdges[best][0]], b = c[kEdges[best][1]];
      if (a > b) std::swap(a, b);
      if (seen.insert({a, b}).second) candidates.push_back({a, b});
    }

    if (candidates.empty()) break;

    SplitResult sr = split_edges(m, candidates);
    res.n_split += sr.n_split;

    // Compose cell_parent across rounds: parent_total[i] =
    // parent_total_prev[result.cell_parent[i]].
    std::vector<int> composed(sr.cell_parent.size());
    for (std::size_t i = 0; i < sr.cell_parent.size(); ++i)
      composed[i] = res.cell_parent[sr.cell_parent[i]];
    res.cell_parent = std::move(composed);

    res.n_smoothed += laplacian_smooth(m, smooth_iters, 0.5);
    m.finalize();

    qs = mesh_quality(m, q_thresh);
  }

  res.min_q_after = mesh_quality(m, q_thresh).min_q;
  return res;
}

}  // namespace cp
