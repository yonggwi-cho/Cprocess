#pragma once
#include <cstdint>
#include <unordered_map>
#include <vector>

#include "mesh.hpp"

namespace cp {

// Connectivity helper for local remeshing and mesh smoothing. Rebuilt from a
// Mesh whenever its topology changes. All ids index into Mesh::nodes / ::cells.
struct MeshTopology {
  // sorted edge (a<b) packed into 64 bits -> incident cell ids
  std::unordered_map<std::uint64_t, std::vector<int>> edge_cells;
  std::vector<std::vector<int>> node_cells;  // node -> incident cell ids
  std::vector<std::vector<int>> node_adj;    // node -> edge-adjacent node ids
  std::vector<char> node_boundary;           // node lies on a boundary face
  std::vector<char> node_interface;          // node touches >1 region tag

  void build(const Mesh& m);

  // Packs an unordered node pair into a key with min in the high word.
  static std::uint64_t edge_key(int a, int b) {
    const std::uint32_t lo = static_cast<std::uint32_t>(a < b ? a : b);
    const std::uint32_t hi = static_cast<std::uint32_t>(a < b ? b : a);
    return (static_cast<std::uint64_t>(lo) << 32) | hi;
  }

  // Cells incident to edge (a,b); empty if the edge is absent.
  const std::vector<int>& edge_incident(int a, int b) const;
};

}  // namespace cp
