#pragma once
#include "mesh.hpp"
#include "topology.hpp"
#include <functional>
#include <vector>

namespace cp {

// Area-weighted average outward normal for each boundary/interface node.
// Interior nodes receive {0,0,0}.
std::vector<Vec3> compute_node_normals(const Mesh& m, const MeshTopology& topo);

struct AleResult {
  int n_moved = 0;     // boundary nodes displaced
  int n_skipped = 0;   // nodes skipped due to inversion guard
};

// Move boundary nodes along their outward normal by v_n*dt.
// mask_fn(node_idx) -> true means the node is masked (blocked from moving).
// Interior nodes are then Laplacian-smoothed (smooth_iters sweeps) to propagate
// the displacement inward without inversions.
// Caller must call Mesh::finalize() afterward to rebuild geometry.
AleResult ale_move(Mesh& m, const MeshTopology& topo,
                   const std::vector<Vec3>& node_normals,
                   double v_n, double dt,
                   std::function<bool(int)> mask_fn = nullptr,
                   int smooth_iters = 3);

// Rescale per-cell concentrations to conserve mass (C*V) when cell volumes
// change due to ALE node movement (topology-invariant).
void rescale_fields_for_volume_change(
    std::vector<std::vector<double>>& fields,
    const std::vector<double>& old_vol,
    const std::vector<double>& new_vol);

}  // namespace cp
