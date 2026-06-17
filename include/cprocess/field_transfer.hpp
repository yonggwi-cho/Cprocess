#pragma once
#include "mesh.hpp"
#include <vector>

namespace cp {

// Transfer a per-cell scalar field from mesh_old to mesh_new using
// nearest-centroid assignment. For each cell in mesh_new the closest
// cell centroid in mesh_old is found and its value is assigned.
std::vector<double> transfer_field_nearest(
    const Mesh& mesh_old, const std::vector<double>& field_old,
    const Mesh& mesh_new);

// Returns relative mass-conservation error:
//   |sum(C_new*V_new) - sum(C_old*V_old)| / sum(C_old*V_old)
double check_mass_conservation(
    const Mesh& mesh_old, const std::vector<double>& field_old,
    const Mesh& mesh_new, const std::vector<double>& field_new);

}  // namespace cp
