#pragma once
#include <vector>

#include "mesh.hpp"

namespace cp {

// Point-in-cell lookup for tetrahedral meshes via a uniform background
// grid of cell bounding boxes. Thread-safe after construction (all queries
// are const), which makes it shareable across Monte Carlo worker threads.
class CellLocator {
 public:
  explicit CellLocator(const Mesh& mesh);

  // Returns the index of a cell containing p (points on shared faces may
  // match either neighbour), or -1 if p lies outside the mesh.
  int locate(const Vec3& p) const;

 private:
  const Mesh& mesh_;
  Vec3 lo_, hi_;
  int nx_ = 1, ny_ = 1, nz_ = 1;
  double ix_ = 0, iy_ = 0, iz_ = 0;  // inverse voxel edge lengths
  std::vector<int> start_;           // voxel -> range start into items_
  std::vector<int> items_;           // candidate cell ids
};

}  // namespace cp
