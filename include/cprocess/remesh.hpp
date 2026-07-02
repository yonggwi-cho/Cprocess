#pragma once
#include <utility>
#include <vector>

#include "mesh.hpp"
#include "vec3.hpp"

namespace cp {

// Normalized tetrahedron quality in (0, 1]: q = 6*sqrt(2)*V / l_rms^3, where
// V is the volume and l_rms the root-mean-square edge length. A regular tet
// scores 1; slivers and degenerate cells approach 0. Sign follows the volume.
double tet_quality(const Vec3& a, const Vec3& b, const Vec3& c, const Vec3& d);

struct QualityStats {
  double min_q = 1.0;   // worst cell
  double mean_q = 0.0;  // average over cells
  int n_sliver = 0;     // cells below sliver_thresh
  int worst_cell = -1;
};

QualityStats mesh_quality(const Mesh& m, double sliver_thresh = 0.1);

// Splits edge (a,b) at its midpoint, subdividing every incident tet into two.
// New cells inherit the region of their parent. Returns the new node id, or
// -1 if the edge is not present. Calls Mesh::finalize() before returning
// (delegates to split_edges).
int edge_split(Mesh& m, int a, int b);

struct SplitResult {
  int n_split = 0;               // number of edges actually split
  int n_skipped = 0;              // edges skipped due to conflicts
  std::vector<int> cell_parent;   // new mesh cell -> old cell index
  int last_new_node = -1;         // node id of the last successfully split edge's midpoint, -1 if none
};

// Splits every edge in `edges` (node id pairs) at its midpoint. Edges whose
// incident cells overlap an already-split edge's incident cells (original or
// new) within this call are skipped (counted in n_skipped) rather than
// applied; callers wanting a full split should retry skipped edges in a
// subsequent pass. Calls Mesh::finalize() once before returning.
SplitResult split_edges(Mesh& m, const std::vector<std::pair<int, int>>& edges);

// Redistributes a per-cell intensive field (e.g. concentration, cm^-3) from
// the old mesh to the new mesh produced by split_edges: each new cell copies
// its parent's value. Intensive quantities are preserved under this copy
// because a split cell's children sum to the parent's volume.
std::vector<double> redistribute_field(const std::vector<double>& conc,
                                       const std::vector<int>& cell_parent);

// Laplacian smoothing of interior nodes; boundary and interface nodes are
// pinned. A candidate move is rejected if it would invert any incident tet.
// Performs `iters` sweeps with relaxation `omega`. Returns the number of node
// moves applied. The caller must call Mesh::finalize() afterward.
int laplacian_smooth(Mesh& m, int iters = 5, double omega = 0.5);

}  // namespace cp
