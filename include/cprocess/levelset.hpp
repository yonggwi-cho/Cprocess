#pragma once
#include "mesh.hpp"
#include "topology.hpp"
#include <vector>

namespace cp {

// Initialize φ as signed distance to interface.
// Cells with region_tag are "inside" (φ < 0), others "outside" (φ > 0).
// φ is defined on nodes (vertex-centered).
std::vector<double> levelset_init(const Mesh& m, int region_tag);

// Same as above but generalized to an arbitrary per-cell inside/outside
// mask (P2-5): inside[ci] != 0 marks cell ci as "inside" (φ < 0). Needed
// because "gas" in this codebase can span several distinct region tags
// (each etch()/deposit() call may mint its own synthetic gas tag), so a
// single region_tag can't identify "everything that isn't solid".
std::vector<double> levelset_init(const Mesh& m, const std::vector<char>& inside);

// First-order upwind advection of φ on tet mesh for one time step.
// v_n: normal velocity at each node (positive = expand interface outward).
// dt: time step (same units as mesh length).
void levelset_advect(const Mesh& m, std::vector<double>& phi,
                     const std::vector<double>& v_n, double dt);

// Node-upwind Godunov advection with CFL substepping (P2-5). Improves on
// the max-|finite-difference| gradient above with the standard first-order
// Hamilton-Jacobi upwind scheme:
//   F_i > 0: |∇φ|_i = max_j max( (φ_i − φ_j)/d_ij, 0 )
//   F_i < 0: |∇φ|_i = max_j max( (φ_j − φ_i)/d_ij, 0 )
// (j ranges over topo.node_adj[i], d_ij the edge length). `dt` is the total
// time to advance; internally sub-stepped at dt_sub = 0.5*h_min/max|F_node|
// (h_min = shortest mesh edge) for stability, so callers don't need to
// worry about CFL themselves.
void levelset_advect(const Mesh& m, const MeshTopology& topo,
                     std::vector<double>& phi,
                     const std::vector<double>& F_node, double dt);

// Reinitialize φ to signed distance via iterative PDE relaxation (Sussman).
// Runs max_iters pseudo-time steps with CFL ~ 0.5.
void levelset_reinit(const Mesh& m, std::vector<double>& phi, int max_iters = 10);

// Update cell_region based on average φ of the cell's nodes.
// Cells with avg φ < 0 → inside_tag, else → outside_tag.
void levelset_update_regions(Mesh& m, const std::vector<double>& phi,
                              int inside_tag, int outside_tag);

}  // namespace cp
