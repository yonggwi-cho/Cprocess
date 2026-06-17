#include "cprocess/levelset.hpp"
#include <algorithm>
#include <cmath>
#include <limits>
#include <queue>

namespace cp {

std::vector<double> levelset_init(const Mesh& m, int region_tag) {
  const int nn = static_cast<int>(m.nodes.size());
  const int nc = static_cast<int>(m.cells.size());

  // Mark nodes by region membership (any incident cell with region_tag → inside).
  std::vector<int> node_inside(nn, 0), node_outside(nn, 0);
  for (int ci = 0; ci < nc; ++ci) {
    const bool inside = (m.cell_region[ci] == region_tag);
    for (int v : m.cells[ci]) {
      if (inside) ++node_inside[v]; else ++node_outside[v];
    }
  }

  // Assign φ = ±1 initially; interface nodes get small magnitude.
  std::vector<double> phi(nn);
  for (int i = 0; i < nn; ++i) {
    if (node_inside[i] > 0 && node_outside[i] > 0) {
      phi[i] = 0.0;  // interface node
    } else if (node_inside[i] > 0) {
      phi[i] = -1.0;
    } else {
      phi[i] = 1.0;
    }
  }

  // BFS to build approximate signed distance from interface.
  // Propagate distance through node adjacency.
  // Build node adjacency on the fly from cells.
  std::vector<std::vector<int>> adj(nn);
  for (int ci = 0; ci < nc; ++ci) {
    const auto& c = m.cells[ci];
    constexpr int edges[6][2] = {{0,1},{0,2},{0,3},{1,2},{1,3},{2,3}};
    for (auto& e : edges) {
      adj[c[e[0]]].push_back(c[e[1]]);
      adj[c[e[1]]].push_back(c[e[0]]);
    }
  }
  // Deduplicate.
  for (int i = 0; i < nn; ++i) {
    std::sort(adj[i].begin(), adj[i].end());
    adj[i].erase(std::unique(adj[i].begin(), adj[i].end()), adj[i].end());
  }

  // Dijkstra-like propagation of distances from interface nodes.
  // Distance stored in phi magnitude; sign from initial assignment.
  std::vector<double> dist(nn, std::numeric_limits<double>::max());
  // Use priority queue: (dist, node).
  using P = std::pair<double, int>;
  std::priority_queue<P, std::vector<P>, std::greater<P>> pq;

  // Seed: interface nodes (phi == 0) or near-interface.
  for (int i = 0; i < nn; ++i) {
    if (phi[i] == 0.0) { dist[i] = 0.0; pq.push({0.0, i}); }
  }

  while (!pq.empty()) {
    auto [d, u] = pq.top(); pq.pop();
    if (d > dist[u]) continue;
    for (int v : adj[u]) {
      const Vec3 dv = m.nodes[v] - m.nodes[u];
      const double edge_len = std::sqrt(dot(dv, dv));
      const double nd = dist[u] + edge_len;
      if (nd < dist[v]) {
        dist[v] = nd;
        pq.push({nd, v});
      }
    }
  }

  // Apply sign and distance.
  for (int i = 0; i < nn; ++i) {
    const double s = (phi[i] <= 0.0) ? -1.0 : 1.0;
    phi[i] = s * (dist[i] < std::numeric_limits<double>::max() ? dist[i] : 1.0);
  }

  return phi;
}

void levelset_advect(const Mesh& m, std::vector<double>& phi,
                     const std::vector<double>& v_n, double dt) {
  const int nn = static_cast<int>(m.nodes.size());

  // Build node adjacency.
  std::vector<std::vector<int>> adj(nn);
  for (const auto& c : m.cells) {
    constexpr int edges[6][2] = {{0,1},{0,2},{0,3},{1,2},{1,3},{2,3}};
    for (auto& e : edges) {
      adj[c[e[0]]].push_back(c[e[1]]);
      adj[c[e[1]]].push_back(c[e[0]]);
    }
  }

  // Upwind: φ_new[i] = φ[i] - dt * |grad φ| * v_n[i]
  // Approximate |grad φ| at node i by max finite difference to neighbors.
  std::vector<double> phi_new(nn);
  for (int i = 0; i < nn; ++i) {
    double grad_mag = 0.0;
    for (int j : adj[i]) {
      const Vec3 dv = m.nodes[j] - m.nodes[i];
      const double len = std::sqrt(dot(dv, dv));
      if (len < 1e-30) continue;
      const double dphi = phi[j] - phi[i];
      // Upwind: use derivative in direction that information comes from.
      const double dphidl = dphi / len;
      grad_mag = std::max(grad_mag, std::fabs(dphidl));
    }
    // Hamilton-Jacobi upwind: φ_t + v_n |∇φ| = 0
    phi_new[i] = phi[i] - dt * v_n[i] * grad_mag;
  }
  phi = std::move(phi_new);
}

void levelset_reinit(const Mesh& m, std::vector<double>& phi, int max_iters) {
  const int nn = static_cast<int>(m.nodes.size());

  // Sussman PDE: φ_τ = sign(φ₀)(1 - |∇φ|)
  // Discretize with a simple node-based upwind.
  std::vector<std::vector<int>> adj(nn);
  std::vector<double> h(nn, std::numeric_limits<double>::max());
  for (const auto& c : m.cells) {
    constexpr int edges[6][2] = {{0,1},{0,2},{0,3},{1,2},{1,3},{2,3}};
    for (auto& e : edges) {
      const int a = c[e[0]], b = c[e[1]];
      adj[a].push_back(b); adj[b].push_back(a);
      const Vec3 dv = m.nodes[b] - m.nodes[a];
      const double len = std::sqrt(dot(dv, dv));
      h[a] = std::min(h[a], len);
      h[b] = std::min(h[b], len);
    }
  }

  const std::vector<double> phi0 = phi;
  for (int it = 0; it < max_iters; ++it) {
    std::vector<double> phi_new(nn);
    for (int i = 0; i < nn; ++i) {
      double grad_mag = 0.0;
      for (int j : adj[i]) {
        const Vec3 dv = m.nodes[j] - m.nodes[i];
        const double len = std::sqrt(dot(dv, dv));
        if (len < 1e-30) continue;
        grad_mag = std::max(grad_mag, std::fabs(phi[j] - phi[i]) / len);
      }
      const double sgn = (phi0[i] >= 0) ? 1.0 : -1.0;
      const double dtau = 0.5 * (h[i] > 0 ? h[i] : 1e-6);
      phi_new[i] = phi[i] + dtau * sgn * (1.0 - grad_mag);
    }
    phi = std::move(phi_new);
  }
}

void levelset_update_regions(Mesh& m, const std::vector<double>& phi,
                              int inside_tag, int outside_tag) {
  const int nc = static_cast<int>(m.cells.size());
  for (int ci = 0; ci < nc; ++ci) {
    double avg = 0.0;
    for (int v : m.cells[ci]) avg += phi[v];
    avg /= 4.0;
    m.cell_region[ci] = (avg < 0.0) ? inside_tag : outside_tag;
  }
}

}  // namespace cp
