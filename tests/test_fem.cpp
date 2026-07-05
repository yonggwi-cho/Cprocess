#include <cstdio>
#include <cmath>

#include "cprocess/fem.hpp"
#include "cprocess/mesh.hpp"
#include "cprocess/mechanics.hpp"
#include "test_util.hpp"

using namespace cp;

namespace {

// von Mises equivalent stress from a fem.hpp-order Voigt stress vector
// (xx,yy,zz,yz,xz,xy).
double von_mises(const double s[6]) {
  const double sxx = s[0], syy = s[1], szz = s[2];
  const double syz = s[3], sxz = s[4], sxy = s[5];
  const double dxx = sxx - syy, dyy = syy - szz, dzz = szz - sxx;
  return std::sqrt(0.5 * (dxx * dxx + dyy * dyy + dzz * dzz) +
                   3.0 * (syz * syz + sxz * sxz + sxy * sxy));
}

}  // namespace

int main() {
  // ---- Test 1: patch test --------------------------------------------------
  // 2x2x2 box mesh, prescribe a linear displacement field u = A*x on every
  // boundary node; the FEM solution on interior nodes must reproduce the
  // same linear field exactly (to solver tolerance), and every element's
  // stress must be uniform (constant-strain elements exactly represent a
  // linear field).
  {
    Mesh m = make_box_mesh(0, 1, 0, 1, 0, 1, 2, 2, 2);
    const int nn = static_cast<int>(m.nodes.size());
    const int nc = static_cast<int>(m.cells.size());

    const double A[3][3] = {{1e-4, 2e-5, 0.0}, {0.0, -3e-5, 1e-5}, {0.0, 0.0, 5e-5}};
    auto u_exact = [&](const Vec3& p) {
      Vec3 u;
      u.x = A[0][0] * p.x + A[0][1] * p.y + A[0][2] * p.z;
      u.y = A[1][0] * p.x + A[1][1] * p.y + A[1][2] * p.z;
      u.z = A[2][0] * p.x + A[2][1] * p.y + A[2][2] * p.z;
      return u;
    };

    std::vector<double> E_cell(nc, 130e3), nu_cell(nc, 0.28);  // E in MPa
    std::vector<double> eps0(nc * 6, 0.0);
    FemProblem prob = fem_assemble(m, E_cell, nu_cell, eps0);

    const BBox bb = m.bbox();
    const double eps = 1e-9 * (bb.hi.x - bb.lo.x);
    for (int i = 0; i < nn; ++i) {
      const Vec3& p = m.nodes[i];
      const bool boundary = std::fabs(p.x - bb.lo.x) < eps || std::fabs(p.x - bb.hi.x) < eps ||
                            std::fabs(p.y - bb.lo.y) < eps || std::fabs(p.y - bb.hi.y) < eps ||
                            std::fabs(p.z - bb.lo.z) < eps || std::fabs(p.z - bb.hi.z) < eps;
      if (!boundary) continue;
      const Vec3 ue = u_exact(p);
      const double comp[3] = {ue.x, ue.y, ue.z};
      for (int d = 0; d < 3; ++d) {
        prob.dof_fixed[3 * i + d] = 1;
        prob.dof_val[3 * i + d] = comp[d];
      }
    }
    fem_apply_bc(prob);

    std::vector<double> u;
    SolveResult res = fem_solve(prob, u, 1e-12, 5000);
    std::printf("patch test: cg converged=%d iters=%d resid=%.3e\n", res.converged,
               res.iters, res.resid);
    CHECK(res.converged);

    double max_rel_err = 0.0, unorm = 0.0;
    for (int i = 0; i < nn; ++i) {
      const Vec3 ue = u_exact(m.nodes[i]);
      const double comp[3] = {ue.x, ue.y, ue.z};
      for (int d = 0; d < 3; ++d) {
        const double err = std::fabs(u[3 * i + d] - comp[d]);
        max_rel_err = std::max(max_rel_err, err);
        unorm = std::max(unorm, std::fabs(comp[d]));
      }
    }
    const double rel = max_rel_err / unorm;
    std::printf("patch test: max |u_fem - u_exact| = %.3e (rel %.3e)\n", max_rel_err, rel);
    CHECK(rel < 1e-10);

    std::vector<double> stress = fem_element_stress(m, u, E_cell, nu_cell, eps0);
    double sxx0 = stress[0];
    double max_diff = 0.0;
    for (int c = 0; c < nc; ++c)
      max_diff = std::max(max_diff, std::fabs(stress[c * 6] - sxx0));
    const double rel_stress = max_diff / std::fabs(sxx0);
    std::printf("patch test: element sigma_xx uniform to rel %.3e\n", rel_stress);
    CHECK(rel_stress < 1e-8);
  }

  std::printf("fem tests passed\n");
  return 0;
}
