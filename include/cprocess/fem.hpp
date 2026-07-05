#pragma once
#include <vector>

#include "mesh.hpp"
#include "sparse.hpp"
#include "vec3.hpp"

namespace cp {

// Linear-elastic FEM on constant-strain (P1) tetrahedra (P2-6).
//
// Voigt convention used throughout this module (note: this differs from the
// (xx,yy,zz,xy,yz,xz) order used by cp::maxwell_update/StressState -- callers
// crossing that boundary must permute components 3..5):
//   index:  0   1   2   3   4   5
//   comp : xx  yy  zz  yz  xz  xy
//
// Units are whatever the caller uses consistently (E and stress in the same
// unit, e.g. both in MPa); lengths follow Mesh's cm convention.

// Shape-function gradients (constant per element) and signed volume of a
// tetrahedron with vertices p[0..3] (already positive-volume-oriented, as
// produced by Mesh::finalize()).
//   grad[0] = (p[3]-p[1]) x (p[2]-p[1]) / (6V)
//   grad[1] = (p[2]-p[0]) x (p[3]-p[0]) / (6V)
//   grad[2] = (p[3]-p[0]) x (p[1]-p[0]) / (6V)
//   grad[3] = (p[1]-p[0]) x (p[2]-p[0]) / (6V)
// Sum grad[0..3] == 0 (asserted internally in debug via the caller's tests).
void tet_shape_grads(const Vec3 p[4], double& vol, Vec3 grad[4]);

// 6x12 strain-displacement matrix from the 4 shape-function gradients.
void tet_B_matrix(const Vec3 grad[4], double B[6][12]);

// 6x6 isotropic elastic stiffness (Voigt, engineering shear strain
// convention matching tet_B_matrix): lambda = E*nu/((1+nu)(1-2nu)),
// mu = E/(2(1+nu)); top-left 3x3 block = lambda + 2*mu*delta_ij, bottom-right
// diagonal (shear rows/cols) = mu, all other entries zero.
void iso_D_matrix(double E, double nu, double D[6][6]);

// Assembled linear system: 3 dof/node (ux,uy,uz), dof id = 3*node + comp,
// expanded into a plain scalar CSR (cg_ilu0 solves this directly).
struct FemProblem {
  CSR K;                        // 3*n_nodes square, SPD after BC application
  std::vector<double> f;        // load vector, size 3*n_nodes
  std::vector<char> dof_fixed;  // 1 = Dirichlet-constrained, size 3*n_nodes
  std::vector<double> dof_val;  // prescribed value for fixed dofs
};

// Assembles K and f (eigenstrain load only; no BC applied yet) for the given
// mesh. E_cell/nu_cell/eps0_cell are per-cell (nc, nc, nc*6 Voigt eigenstrain
// in this module's index order). dof_fixed/dof_val come back sized and
// zeroed (unconstrained); set entries and call fem_apply_bc() before solving.
FemProblem fem_assemble(const Mesh& m, const std::vector<double>& E_cell,
                        const std::vector<double>& nu_cell,
                        const std::vector<double>& eps0_cell);

// Applies Dirichlet BCs already recorded in p.dof_fixed/p.dof_val via the
// standard diagonal-fixing / symmetric-elimination procedure: for a fixed
// dof i, its row becomes the identity row (rhs = prescribed value); for every
// other row j with a nonzero column i, K(j,i)*val[i] is subtracted from f[j]
// before column i is zeroed (off the diagonal). This preserves SPD-ness of
// the reduced system for the free dofs while keeping K's global size fixed.
void fem_apply_bc(FemProblem& p);

// Solves p.K u = p.f via cg_ilu0 (p must already have BCs applied). Returns
// the CG convergence diagnostics; u comes back sized 3*n_nodes.
SolveResult fem_solve(const FemProblem& p, std::vector<double>& u,
                      double rtol = 1e-10, int maxit = 5000);

// Per-cell Voigt stress sigma = D*(B*u_local - eps0), same index order as
// eps0_cell above. Returned as nc*6, row-major per cell.
std::vector<double> fem_element_stress(const Mesh& m,
                                       const std::vector<double>& u,
                                       const std::vector<double>& E_cell,
                                       const std::vector<double>& nu_cell,
                                       const std::vector<double>& eps0_cell);

}  // namespace cp
