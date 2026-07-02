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

struct RepairResult {
  int n_smoothed = 0;     // nodes moved by smoothing (summed over all passes)
  int n_split = 0;        // edges split (summed over all rounds)
  int n_flips = 0;        // 2-3 + 3-2 flips applied (summed over all rounds)
  double min_q_before = 0, min_q_after = 0;
  std::vector<int> cell_parent;  // maps final mesh cell -> original cell index
                                  // (identity if no splits occurred). CLEARED
                                  // (empty) when flips occurred — merged cells
                                  // have no single parent; pass fields to
                                  // repair_quality for correct transfer.
};

// Result of a single flip_repair pass.
struct FlipResult { int n_flip23 = 0, n_flip32 = 0; };

// Describes one flip's effect on the cell array: `old_cells` (pre-pass
// indices) were consumed and replaced by `new_cells` (post-pass, post-
// compaction indices).
struct FlipRemap { std::vector<int> old_cells, new_cells; };

// Attempts local 2-3 (face) and 3-2 (edge) flips on cells with quality
// q < q_thresh, worst first. Both flips preserve the union polyhedron, so the
// boundary surface and region interfaces are never altered: a 2-3 flip only
// targets a face shared by two same-region cells (an interior face), a 3-2
// flip only targets an edge whose incident-cell ring is closed (exactly three
// same-region cells, the ring's other vertices forming one triangle), and
// every flip must conserve group volume to 1e-9 relative (this rejects
// non-convex configurations and any boundary-altering move). A flip is only
// applied if it strictly improves the minimum quality among the cells it
// touches. Runs a single pass (topology is built
// once at entry) and calls Mesh::finalize() before returning. If `remaps` is
// non-null, one FlipRemap is appended per successful flip. If `fields` is
// non-null, every listed field vector is resized/rewritten in place so that
// each final cell's value is either copied from its untouched source cell or
// set to the volume-weighted average (using pre-flip cell_vol) of its flip
// group's old cells.
FlipResult flip_repair(Mesh& m, double q_thresh,
                       std::vector<FlipRemap>* remaps = nullptr,
                       std::vector<std::vector<double>*>* fields = nullptr);

// Repairs cells with quality q < q_thresh by alternating Laplacian smoothing,
// longest-edge splitting and local 2-3/3-2 flipping, up to `max_rounds`
// rounds. Each round: smooth `smooth_iters` sweeps, recheck quality (stop
// early if no slivers remain), split the longest edge of every remaining
// sliver cell, run flip_repair, and smooth again. Never throws; results
// (including any residual slivers) are reported via the return value. Does
// not call ale_move; the caller is responsible for composing this with any
// prior mesh motion. If `fields` is non-null, every listed per-cell field is
// carried through splits (copy from parent) and flips (volume-weighted
// average) in place.
RepairResult repair_quality(Mesh& m, std::vector<std::vector<double>*>* fields,
                            double q_thresh = 0.1, int max_rounds = 3,
                            int smooth_iters = 5);

// Backward-compatible overload without field transfer.
RepairResult repair_quality(Mesh& m, double q_thresh = 0.1, int max_rounds = 3,
                            int smooth_iters = 5);

}  // namespace cp
