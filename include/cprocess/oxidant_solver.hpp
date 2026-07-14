#pragma once
#include <iosfwd>
#include <utility>
#include <vector>

#include "mesh.hpp"

namespace cp {

// Result of a steady-state oxidant diffusion-reaction solve (P2-4).
struct OxidantResult {
  std::vector<double> conc;  // per-cell oxidant concentration [cm^-3]; 0 outside the active domain
  // Per-interface-face growth velocity: (face index, v_n [cm/s], Si-ward positive).
  std::vector<std::pair<int, double>> iface_vn;
};

// Steady-state div(D grad C) = 0 over the "active" cell set (typically the
// oxide sub-mesh; may also include nitride cells with a near-zero D so a
// small amount of oxidant can bleed laterally under a mask edge -- see the
// two-argument-D overload below).
//
//   gas/oxide surface (an active cell's boundary face on patch "zmax"):
//     Dirichlet C = c_gas
//   active/Si interface face (active cell adjacent to a `si_mask` cell):
//     Robin: reaction flux F = ks * C_face, C_face eliminated via the
//     series resistance t_robin = 1 / (1/(d_ox*g) + 1/(ks*area)); this
//     acts as a sink added to the cell's diagonal (F_f = t_robin * C_owner).
//   all other faces (side walls, active/active-inactive-non-Si, nitride's
//     own top/side boundaries): zero flux.
//
// `d_ox` is a uniform oxidant diffusivity for every active cell. Cells with
// active_mask == 0 get an identity row (C == 0) -- consistent with "an
// empty/inactive cell has no oxidant".  SPD system, solved with cg_ilu0.
//
// Implementation note: `m` must be a make_box_mesh-style Kuhn-tet mesh (6
// tets per hex, emitted in hex-scan order -- true for every mesh this
// codebase builds/extends). The solve is done at hex-aggregate resolution
// (see src/oxidant_solver.cpp) to sidestep the severe non-orthogonality of
// the internal per-hex diagonal tet faces; `conc` still returns one value
// per tet (each tet copies its hex's value).
OxidantResult solve_oxidant(const Mesh& m, const std::vector<char>& oxide_mask,
                            const std::vector<char>& si_mask, double d_ox,
                            double ks, double c_gas, std::ostream* log);

// General form with a per-cell diffusivity (P2-4 Stage B): lets nitride
// cells participate in the active domain at D_nitride << D_ox, which is
// what allows a small amount of lateral oxidant leakage under the mask
// edge (combined with lateral diffusion from the open field) to produce a
// bird's-beak taper. `active_mask` marks every cell in the linear system's
// non-trivial rows (oxide U nitride, typically); `d_cell` gives each such
// cell's diffusivity (0 elsewhere is fine, ignored for inactive cells).
// P3-e: optional per-face reaction-rate multiplier for the Robin (active/Si)
// faces, indexed by face id (size m.faces.size()). nullptr (default) leaves
// every Robin face's t_robin at its pre-P3-e value (ks unscaled) -- bit-
// identical to the 5-argument call above. When non-null, `ks_face_scale[fi]`
// multiplies `ks` for that face only (e.g. stress-slowed interface reaction,
// scale in (0, 1]); faces never visited by the Robin branch ignore their
// entry.
OxidantResult solve_oxidant(const Mesh& m, const std::vector<double>& d_cell,
                            const std::vector<char>& active_mask,
                            const std::vector<char>& si_mask, double ks,
                            double c_gas, std::ostream* log,
                            const std::vector<double>* ks_face_scale = nullptr);

}  // namespace cp
