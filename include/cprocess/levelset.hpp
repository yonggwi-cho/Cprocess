#pragma once
#include "mesh.hpp"
#include <vector>

namespace cp {

// Initialize φ as signed distance to interface.
// Cells with region_tag are "inside" (φ < 0), others "outside" (φ > 0).
// φ is defined on nodes (vertex-centered).
std::vector<double> levelset_init(const Mesh& m, int region_tag);

// First-order upwind advection of φ on tet mesh for one time step.
// v_n: normal velocity at each node (positive = expand interface outward).
// dt: time step (same units as mesh length).
void levelset_advect(const Mesh& m, std::vector<double>& phi,
                     const std::vector<double>& v_n, double dt);

// Reinitialize φ to signed distance via iterative PDE relaxation (Sussman).
// Runs max_iters pseudo-time steps with CFL ~ 0.5.
void levelset_reinit(const Mesh& m, std::vector<double>& phi, int max_iters = 10);

// Update cell_region based on average φ of the cell's nodes.
// Cells with avg φ < 0 → inside_tag, else → outside_tag.
void levelset_update_regions(Mesh& m, const std::vector<double>& phi,
                              int inside_tag, int outside_tag);

}  // namespace cp
