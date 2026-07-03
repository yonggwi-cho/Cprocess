#include "cprocess/remesh.hpp"

#include <algorithm>
#include <cmath>
#include <numeric>
#include <set>
#include <stdexcept>

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

namespace {

// Orients (v0,v1,v2,v3) to positive signed volume (mirroring Mesh::finalize's
// convention) and returns false if the resulting tet is degenerate.
bool make_oriented_tet(const Mesh& m, int v0, int v1, int v2, int v3,
                       std::array<int, 4>& out, double& q, double& vol) {
  std::array<int, 4> t{v0, v1, v2, v3};
  double v = signed_vol(m.nodes[t[0]], m.nodes[t[1]], m.nodes[t[2]], m.nodes[t[3]]);
  if (v < 0) { std::swap(t[2], t[3]); v = -v; }
  if (v <= 1e-30) return false;
  out = t;
  q = tet_quality(m.nodes[t[0]], m.nodes[t[1]], m.nodes[t[2]], m.nodes[t[3]]);
  vol = v;
  return true;
}

}  // namespace

FlipResult flip_repair(Mesh& m, double q_thresh, std::vector<FlipRemap>* remaps,
                       std::vector<std::vector<double>*>* fields) {
  FlipResult res;
  const int nc0 = static_cast<int>(m.cells.size());
  if (nc0 == 0) return res;

  const std::vector<double> vol0 = m.cell_vol;  // snapshot before any mutation

  MeshTopology topo;
  topo.build(m);

  // Pre-pass quality snapshot (node positions never change in this pass, so
  // this stays valid for any untouched cell throughout).
  std::vector<double> qual(nc0);
  for (int ci = 0; ci < nc0; ++ci) {
    const auto& c = m.cells[ci];
    qual[ci] = tet_quality(m.nodes[c[0]], m.nodes[c[1]], m.nodes[c[2]], m.nodes[c[3]]);
  }

  // Face-neighbour list per cell, restricted to internal, same-region faces
  // (region interfaces are never flip candidates).
  struct FaceNb { int nb; std::array<int, 3> tri; };
  std::vector<std::vector<FaceNb>> cell_face_nb(nc0);
  for (const auto& f : m.faces) {
    if (f.neigh < 0) continue;
    if (m.cell_region[f.owner] != m.cell_region[f.neigh]) continue;
    cell_face_nb[f.owner].push_back({f.neigh, f.n});
    cell_face_nb[f.neigh].push_back({f.owner, f.n});
  }

  std::vector<int> order;
  for (int ci = 0; ci < nc0; ++ci) if (qual[ci] < q_thresh) order.push_back(ci);
  std::sort(order.begin(), order.end(), [&](int a, int b) { return qual[a] < qual[b]; });

  std::vector<char> cell_dead(nc0, 0);
  std::vector<int> slot_group(nc0, -1);
  std::vector<std::vector<int>> group_old;    // group id -> old (pre-pass) cell ids
  std::vector<std::vector<int>> group_slots;  // group id -> current slot ids

  auto alive = [&](int ci) {
    return ci >= 0 && static_cast<std::size_t>(ci) < cell_dead.size() &&
           !cell_dead[ci] && slot_group[ci] == -1;
  };

  for (int ci : order) {
    if (!alive(ci)) continue;
    bool applied = false;

    // ---- 3-2 (edge collapse of 3 tets sharing an edge) ---------------------
    for (int e = 0; e < 6 && !applied; ++e) {
      const int a = m.cells[ci][kEdges[e][0]], b = m.cells[ci][kEdges[e][1]];
      const std::vector<int>& inc = topo.edge_incident(a, b);
      if (inc.size() != 3) continue;
      const int c0 = inc[0], c1 = inc[1], c2 = inc[2];
      if (!alive(c0) || !alive(c1) || !alive(c2)) continue;
      if (m.cell_region[c0] != m.cell_region[c1] ||
          m.cell_region[c0] != m.cell_region[c2]) continue;
      // Note: no explicit node_boundary(a)/node_boundary(b) check here — the
      // ring-closure verification below (three tets forming a closed
      // triangular fan around the axis edge) already rejects open boundary
      // fans, which is the only way a mesh-boundary edge could report
      // exactly 3 incident cells.

      auto ring_of = [&](int cc, int& r0, int& r1) {
        int k = 0; int rr[2] = {-1, -1};
        for (int v : m.cells[cc]) if (v != a && v != b) rr[k++] = v;
        r0 = rr[0]; r1 = rr[1];
      };
      int p0, q0, p1, q1, p2, q2;
      ring_of(c0, p0, q0); ring_of(c1, p1, q1); ring_of(c2, p2, q2);

      const int x1 = p0, x2 = q0;
      int x3 = -1, second = -1, third = -1;
      if (p1 == x2 || q1 == x2) { second = c1; x3 = (p1 == x2) ? q1 : p1; third = c2; }
      else if (p2 == x2 || q2 == x2) { second = c2; x3 = (p2 == x2) ? q2 : p2; third = c1; }
      else continue;
      (void)second;
      int tp, tq;
      ring_of(third, tp, tq);
      if (!((tp == x3 && tq == x1) || (tp == x1 && tq == x3))) continue;

      std::array<int, 4> t1, t2;
      double q1v, q2v, v1v, v2v;
      if (!make_oriented_tet(m, x1, x2, x3, a, t1, q1v, v1v)) continue;
      if (!make_oriented_tet(m, x1, x2, x3, b, t2, q2v, v2v)) continue;
      const double newQ = std::min(q1v, q2v);
      const double oldQ = std::min({qual[c0], qual[c1], qual[c2]});
      if (newQ <= oldQ) continue;
      // Volume conservation rejects non-convex (reflex-edge) configurations
      // where positively-oriented replacement tets would overlap.
      const double volOld = vol0[c0] + vol0[c1] + vol0[c2];
      if (std::fabs((v1v + v2v) - volOld) > 1e-9 * volOld) continue;

      const int g = static_cast<int>(group_old.size());
      group_old.push_back({c0, c1, c2});
      group_slots.push_back({});
      m.cells[c0] = t1; slot_group[c0] = g; group_slots[g].push_back(c0);
      m.cells[c1] = t2; slot_group[c1] = g; group_slots[g].push_back(c1);
      cell_dead[c2] = 1; slot_group[c2] = g;

      ++res.n_flip32;
      applied = true;
    }
    if (applied) continue;

    // ---- 2-3 (face split into 3 tets) --------------------------------------
    for (const auto& fnb : cell_face_nb[ci]) {
      const int nb = fnb.nb;
      if (!alive(ci) || !alive(nb)) continue;
      const int a = fnb.tri[0], b = fnb.tri[1], c = fnb.tri[2];
      // Note: cell_face_nb only contains internal (neigh >= 0), same-region
      // faces, so this candidate face is never a mesh-boundary or interface
      // face; the exterior surface is unaffected regardless of whether a, b,
      // c also happen to sit on the boundary elsewhere.

      int d = -1, e2 = -1;
      for (int v : m.cells[ci]) if (v != a && v != b && v != c) d = v;
      for (int v : m.cells[nb]) if (v != a && v != b && v != c) e2 = v;
      if (d < 0 || e2 < 0) continue;

      std::array<int, 4> t1, t2, t3;
      double q1v, q2v, q3v, v1v, v2v, v3v;
      if (!make_oriented_tet(m, a, b, d, e2, t1, q1v, v1v)) continue;
      if (!make_oriented_tet(m, b, c, d, e2, t2, q2v, v2v)) continue;
      if (!make_oriented_tet(m, c, a, d, e2, t3, q3v, v3v)) continue;
      const double newQ = std::min({q1v, q2v, q3v});
      const double oldQ = std::min(qual[ci], qual[nb]);
      if (newQ <= oldQ) continue;
      // Volume conservation rejects non-convex configurations where
      // positively-oriented replacement tets would overlap.
      const double volOld = vol0[ci] + vol0[nb];
      if (std::fabs((v1v + v2v + v3v) - volOld) > 1e-9 * volOld) continue;

      const int g = static_cast<int>(group_old.size());
      group_old.push_back({ci, nb});
      group_slots.push_back({});
      m.cells[ci] = t1; slot_group[ci] = g; group_slots[g].push_back(ci);
      m.cells[nb] = t2; slot_group[nb] = g; group_slots[g].push_back(nb);
      const int newIdx = static_cast<int>(m.cells.size());
      m.cells.push_back(t3);
      m.cell_region.push_back(m.cell_region[ci]);
      cell_dead.push_back(0);
      slot_group.push_back(g);
      group_slots[g].push_back(newIdx);

      ++res.n_flip23;
      applied = true;
      break;
    }
  }

  // Compact: drop dead slots, keep relative order.
  const int nraw = static_cast<int>(m.cells.size());
  std::vector<int> new_index(nraw, -1);
  std::vector<std::array<int, 4>> new_cells;
  std::vector<int> new_region;
  new_cells.reserve(nraw);
  new_region.reserve(nraw);
  for (int s = 0; s < nraw; ++s) {
    if (cell_dead[s]) continue;
    new_index[s] = static_cast<int>(new_cells.size());
    new_cells.push_back(m.cells[s]);
    new_region.push_back(m.cell_region[s]);
  }

  if (fields) {
    const int newN = static_cast<int>(new_cells.size());
    for (auto* vecp : *fields) {
      if (!vecp) continue;
      std::vector<double>& vec = *vecp;
      std::vector<double> out(newN, 0.0);
      std::vector<double> group_avg(group_old.size(), 0.0);
      std::vector<char> group_ready(group_old.size(), 0);
      for (int s = 0; s < nraw; ++s) {
        const int nidx = new_index[s];
        if (nidx < 0) continue;
        if (slot_group[s] == -1) {
          out[nidx] = (s < static_cast<int>(vec.size())) ? vec[s] : 0.0;
        } else {
          const int g = slot_group[s];
          if (!group_ready[g]) {
            double num = 0, den = 0;
            for (int oc : group_old[g]) { num += vec[oc] * vol0[oc]; den += vol0[oc]; }
            group_avg[g] = den > 0 ? num / den : 0.0;
            group_ready[g] = 1;
          }
          out[nidx] = group_avg[g];
        }
      }
      vec = std::move(out);
    }
  }

  if (remaps) {
    for (std::size_t g = 0; g < group_old.size(); ++g) {
      FlipRemap fr;
      fr.old_cells = group_old[g];
      for (int s : group_slots[g]) fr.new_cells.push_back(new_index[s]);
      remaps->push_back(std::move(fr));
    }
  }

  m.cells = std::move(new_cells);
  m.cell_region = std::move(new_region);
  m.finalize();
  return res;
}

RepairResult repair_quality(Mesh& m, std::vector<std::vector<double>*>* fields,
                            double q_thresh, int max_rounds, int smooth_iters) {
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

    if (!candidates.empty()) {
      SplitResult sr = split_edges(m, candidates);
      res.n_split += sr.n_split;

      // Compose cell_parent across rounds: parent_total[i] =
      // parent_total_prev[result.cell_parent[i]]. Once cell_parent has been
      // cleared by an earlier flip, it stays cleared (no single parent).
      if (!res.cell_parent.empty()) {
        std::vector<int> composed(sr.cell_parent.size());
        for (std::size_t i = 0; i < sr.cell_parent.size(); ++i)
          composed[i] = res.cell_parent[sr.cell_parent[i]];
        res.cell_parent = std::move(composed);
      }

      if (fields) {
        for (auto* vecp : *fields) {
          if (!vecp) continue;
          *vecp = redistribute_field(*vecp, sr.cell_parent);
        }
      }
    }

    // Local 2-3/3-2 flips (never touching interfaces or boundary). Flips
    // merge/split cells, so a single-parent mapping is no longer meaningful:
    // per RepairResult::cell_parent's contract, clear it once any flip
    // fires (field transfer through flips is handled directly by
    // flip_repair via `fields`).
    {
      FlipResult fr = flip_repair(m, q_thresh, nullptr, fields);
      res.n_flips += fr.n_flip23 + fr.n_flip32;
      if (fr.n_flip23 + fr.n_flip32 > 0) res.cell_parent.clear();
    }

    res.n_smoothed += laplacian_smooth(m, smooth_iters, 0.5);
    m.finalize();

    qs = mesh_quality(m, q_thresh);
  }

  res.min_q_after = mesh_quality(m, q_thresh).min_q;
  return res;
}

RepairResult repair_quality(Mesh& m, double q_thresh, int max_rounds,
                            int smooth_iters) {
  return repair_quality(m, nullptr, q_thresh, max_rounds, smooth_iters);
}

RefineResult refine_gradient(Mesh& m, std::vector<std::vector<double>*>& fields,
                             int key_index, double rel_grad_thresh,
                             int max_passes, double max_growth) {
  if (key_index < 0 || static_cast<std::size_t>(key_index) >= fields.size())
    throw std::invalid_argument("refine_gradient: key_index out of range");

  RefineResult res;
  res.n_cells_before = static_cast<int>(m.cells.size());
  const int max_cells = static_cast<int>(res.n_cells_before * max_growth);

  for (int pass = 0; pass < max_passes; ++pass) {
    if (static_cast<int>(m.cells.size()) >= max_cells) break;
    const std::vector<double>& conc = *fields[key_index];
    double global_max = 0.0;
    for (double v : conc) if (v > global_max) global_max = v;
    if (global_max <= 0) break;
    const double floor_val = 1e-3 * global_max;

    std::vector<std::pair<int, int>> candidates;
    std::set<std::pair<int, int>> seen;
    for (const auto& f : m.faces) {
      if (f.neigh < 0) continue;
      const double c_o = conc[f.owner], c_n = conc[f.neigh];
      const double thresh = rel_grad_thresh * std::max({c_o, c_n, floor_val});
      if (std::fabs(c_o - c_n) <= thresh) continue;

      const auto& c = m.cells[f.owner];
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

    // Apply the batch one edge at a time so a single pass cannot blow past
    // max_cells (a face-driven candidate set can be a large fraction of all
    // edges when rel_grad_thresh is small). Before each split, use a fresh
    // topology to know exactly how many cells that split will add, and skip
    // (not abort) edges that would push the total at or beyond the cap --
    // other, cheaper edges later in the batch may still fit.
    for (const auto& edge : candidates) {
      if (static_cast<int>(m.cells.size()) >= max_cells) break;

      MeshTopology etopo;
      etopo.build(m);
      const std::size_t inc = etopo.edge_incident(edge.first, edge.second).size();
      if (inc == 0) continue;
      if (static_cast<int>(m.cells.size() + inc) > max_cells) continue;

      SplitResult sr = split_edges(m, {edge});
      res.n_split_total += sr.n_split;

      for (auto* vecp : fields) {
        if (!vecp) continue;
        *vecp = redistribute_field(*vecp, sr.cell_parent);
      }
    }
    ++res.n_passes;

    if (static_cast<int>(m.cells.size()) >= max_cells) break;
  }

  res.n_cells_after = static_cast<int>(m.cells.size());
  return res;
}

}  // namespace cp
