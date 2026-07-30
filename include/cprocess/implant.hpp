#pragma once
#include <vector>

#include "materials.hpp"
#include "mesh.hpp"

namespace cp {

// Analytic vertical Gaussian implant through the top surface (+z plane of
// the mesh bounding box; a planar top surface is assumed).
//
//   C(d) = dose / (sqrt(2*pi) * dRp) * exp(-(d - Rp)^2 / (2 dRp^2)),
//   d = z_top - z.
//
// With a mask window [x1,x2]x[y1,y2] the lateral spread uses the standard
// error-function convolution with lateral straggle dRl.
struct ImplantParams {
  enum class Profile { gauss, pearson4 };
  const Dopant* dopant = nullptr;
  double dose = 0;          // cm^-2
  double rp = 0, drp = 0;   // cm
  double drl = 0;           // cm; <=0 selects 0.8*dRp
  bool has_window = false;
  double x1 = 0, x2 = 0, y1 = 0, y2 = 0;
  Profile profile = Profile::gauss;
  double gamma = 0.0, beta = 3.0;  // Pearson-IV moments (profile=pearson4)
};

// Adds the implant profile to `conc` for cells where mask is true.
// Returns the number of implanted atoms (integral of the added profile).
//
// W-8: `depth_shift`, when non-null, is a per-cell correction added to the
// depth argument d = z_top - z before the profile is evaluated. It encodes
// column-wise screening by overlying non-silicon layers: for each layer of
// thickness t above the cell the shift contributes (S_layer/S_Si - 1)*t
// (Si-equivalent thickness minus the physical thickness already contained
// in d). Cells whose profile falls inside the screen layer simply receive
// the deep (attenuated) part of the profile — the screened-off dose is NOT
// renormalized into silicon (physical dose loss, matching the MC behavior).
// Passing nullptr reproduces the legacy behavior bit-identically.
double apply_implant(const Mesh& mesh, const std::vector<char>& mask,
                     const ImplantParams& p, std::vector<double>& conc,
                     const std::vector<double>* depth_shift = nullptr);

}  // namespace cp
