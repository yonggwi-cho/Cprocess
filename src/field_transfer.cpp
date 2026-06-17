#include "cprocess/field_transfer.hpp"
#include <cmath>
#include <limits>

namespace cp {

std::vector<double> transfer_field_nearest(
    const Mesh& mesh_old, const std::vector<double>& field_old,
    const Mesh& mesh_new) {
  const int nc_new = static_cast<int>(mesh_new.cells.size());
  const int nc_old = static_cast<int>(mesh_old.cells.size());
  std::vector<double> result(nc_new, 0.0);

  for (int i = 0; i < nc_new; ++i) {
    const Vec3& c = mesh_new.cell_cent[i];
    double best_d2 = std::numeric_limits<double>::max();
    int best_j = 0;
    for (int j = 0; j < nc_old; ++j) {
      const Vec3 d = c - mesh_old.cell_cent[j];
      const double d2 = dot(d, d);
      if (d2 < best_d2) { best_d2 = d2; best_j = j; }
    }
    result[i] = field_old[best_j];
  }
  return result;
}

double check_mass_conservation(
    const Mesh& mesh_old, const std::vector<double>& field_old,
    const Mesh& mesh_new, const std::vector<double>& field_new) {
  double mass_old = 0, mass_new = 0;
  for (std::size_t i = 0; i < field_old.size(); ++i)
    mass_old += field_old[i] * mesh_old.cell_vol[i];
  for (std::size_t i = 0; i < field_new.size(); ++i)
    mass_new += field_new[i] * mesh_new.cell_vol[i];
  if (mass_old <= 0) return 0.0;
  return std::fabs(mass_new - mass_old) / mass_old;
}

}  // namespace cp
