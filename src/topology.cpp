#include "cprocess/topology.hpp"

#include <algorithm>

namespace cp {

namespace {
// The six edges of a tet as local-vertex index pairs.
constexpr int kEdges[6][2] = {{0, 1}, {0, 2}, {0, 3}, {1, 2}, {1, 3}, {2, 3}};
}  // namespace

void MeshTopology::build(const Mesh& m) {
  const int nn = static_cast<int>(m.nodes.size());
  const int nc = static_cast<int>(m.cells.size());

  edge_cells.clear();
  edge_cells.reserve(static_cast<std::size_t>(nc) * 6);
  node_cells.assign(nn, {});
  node_adj.assign(nn, {});
  node_boundary.assign(nn, 0);
  node_interface.assign(nn, 0);

  for (int ci = 0; ci < nc; ++ci) {
    const auto& c = m.cells[ci];
    for (int v = 0; v < 4; ++v) node_cells[c[v]].push_back(ci);
    for (const auto& e : kEdges)
      edge_cells[edge_key(c[e[0]], c[e[1]])].push_back(ci);
  }

  // Node adjacency from unique edges.
  for (const auto& [key, cells] : edge_cells) {
    (void)cells;
    const int a = static_cast<int>(key >> 32);
    const int b = static_cast<int>(key & 0xffffffffu);
    node_adj[a].push_back(b);
    node_adj[b].push_back(a);
  }

  // Boundary nodes: any node on a face with no neighbour.
  for (const auto& f : m.faces)
    if (f.neigh < 0)
      for (int v : f.n) node_boundary[v] = 1;

  // Interface nodes: incident cells span more than one region tag.
  for (int i = 0; i < nn; ++i) {
    if (node_cells[i].empty()) continue;
    const int r0 = m.cell_region[node_cells[i].front()];
    for (int ci : node_cells[i])
      if (m.cell_region[ci] != r0) { node_interface[i] = 1; break; }
  }
}

const std::vector<int>& MeshTopology::edge_incident(int a, int b) const {
  static const std::vector<int> kEmpty;
  auto it = edge_cells.find(edge_key(a, b));
  return it == edge_cells.end() ? kEmpty : it->second;
}

}  // namespace cp
