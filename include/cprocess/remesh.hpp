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

struct RefineResult {
  int n_passes = 0;        // number of passes executed
  int n_split_total = 0;   // total edges split
  int n_cells_before = 0, n_cells_after = 0;
};

// Gradient-driven adaptive refinement. `fields` holds every per-cell field
// that must be carried through the splits (redistribute_field, exact parent
// copy); `key_index` selects which entry of `fields` drives the gradient
// indicator (throws std::invalid_argument if out of range). Each pass scans
// internal faces; a face is a "steep" candidate when
//   |C_owner - C_neigh| > rel_grad_thresh * max(C_owner, C_neigh, 1e-3*global_max)
// where global_max is the current max of the key field. For every steep face,
// the owner cell's longest edge is added (deduplicated) to the split batch.
// Stops when no candidates remain, `max_passes` is reached, or the cell count
// would exceed n_cells_before * max_growth.
RefineResult refine_gradient(Mesh& m, std::vector<std::vector<double>*>& fields,
                             int key_index, double rel_grad_thresh,
                             int max_passes, double max_growth = 4.0);

// Direction-aligned (anisotropic) refinement (M-6). Same pass structure as
// refine_gradient (M-2) -- fields carried through by parent-cell copy,
// max_growth caps total cell count -- but the edge indicator is alignment
// with `direction` times edge length rather than a field gradient:
//   score(edge) = |dot(e_hat, d_hat)| * len(edge)
// (e_hat: unit edge vector, d_hat: direction normalized). An edge is a split
// candidate when score > align_thresh * L_ref, where L_ref is the mean
// length of "direction-aligned" edges in the current mesh (|dot(e_hat,d_hat)|
// > 0.9), falling back to the mean length of all edges if none qualify.
// Candidates are additionally required to satisfy |dot(e_hat,d_hat)| > 0.9
// (the same cutoff used to define L_ref's aligned population) before
// being sorted by score descending (longest, most-aligned first) and handed
// to split_edges (conflict skips are that function's responsibility, same as
// M-2). Measured deviation from a literal reading of the spec: score alone
// is algebraically just dot(e, d_hat) (the edge vector's projection onto
// direction), so on a Kuhn-triangulated box mesh every one-layer-tall edge --
// pure axis edges as well as face/main diagonals that also span one cell in
// the aligned direction -- scores identically; without the extra alignment
// gate, splitting would inject as much lateral (off-axis) structure as
// intended depth resolution (measured on a 6x6x6 unit box: 720 diagonal
// candidates vs. 294 true z-aligned ones, scoring equal). Throws
// std::invalid_argument if |direction| <= 0.
//
// Note: a midpoint split halves an edge's length, so after a few passes a
// once-aligned edge's score falls below align_thresh * L_ref (L_ref itself
// shrinks pass to pass as the aligned population gets shorter) and
// refinement self-terminates -- this is the mechanism that thins layers
// along `direction` without a separate stopping heuristic. max_passes is
// therefore the practical depth control.
RefineResult refine_anisotropic(Mesh& m, const Vec3& direction,
                                std::vector<std::vector<double>*>* fields,
                                double align_thresh, int max_passes,
                                double max_growth = 4.0);

struct CoarsenResult {
  int n_collapsed = 0;      // number of edge collapses actually applied
  int n_rejected = 0;       // candidates rejected (protection/quality/inversion)
  int n_cells_removed = 0;  // cells dropped from the mesh
};

// Collapses each edge (a,b) in `edges` by merging node b into node a (a stays
// put). Boundary and region-interface nodes are never touched (an edge with
// either endpoint protected is rejected). A candidate collapse is simulated
// before being applied: cells incident to both a and b vanish (positive
// volume sums to zero across an internal edge collapse); cells containing
// only b are reconnected with b -> a. If any reconnected cell would end up
// with non-positive signed volume or tet_quality < 0.05, the whole collapse
// is rejected (this substitutes for a full topological link condition). If
// `fields` is non-null, every listed per-cell field is transferred so that
// total mass (Sum C*V) is conserved: reconnected cells are rescaled to keep
// their own mass constant (C_new = C_old * V_old / V_new), and each removed
// cell's mass is split evenly across the surviving cells that share at least
// 3 nodes with it post-collapse (falling back to the collapse's reconnected
// cells if none qualify). Calls Mesh::finalize() once at the end. Orphaned
// nodes (b) are left in place, unreferenced; node array compaction is out of
// scope for this function.
CoarsenResult coarsen(Mesh& m, const std::vector<std::pair<int, int>>& edges,
                      std::vector<std::vector<double>*>* fields = nullptr);

// Selects edge-collapse candidates for `coarsen`. Scans every internal face;
// a face is "low-gradient" when
//   |C_owner - C_neigh| < rel_grad_thresh * max(C_owner, C_neigh, 1e-3*global_max)
// where global_max is the max of `conc`. For each low-gradient face, the
// owner cell's shortest edge is a candidate (deduplicated, both endpoints
// required to be non-boundary and non-interface). Candidates are sorted by
// length ascending and greedily accepted, skipping any whose incident cells
// overlap an already-accepted candidate's incident cells (same conflict
// avoidance as split_edges). At most
// `edge_cells.size() * max_fraction` edges are returned.
std::vector<std::pair<int, int>> select_coarsen_edges(
    const Mesh& m, const std::vector<double>& conc, double rel_grad_thresh,
    double max_fraction = 0.1);

}  // namespace cp
