#include "cprocess/fem.hpp"

#include <algorithm>
#include <cassert>
#include <cmath>

#include "cprocess/topology.hpp"

namespace cp {

void tet_shape_grads(const Vec3 p[4], double& vol, Vec3 grad[4]) {
  const Vec3& p0 = p[0];
  const Vec3& p1 = p[1];
  const Vec3& p2 = p[2];
  const Vec3& p3 = p[3];
  vol = dot(p1 - p0, cross(p2 - p0, p3 - p0)) / 6.0;
  const double inv6v = 1.0 / (6.0 * vol);
  grad[0] = inv6v * cross(p3 - p1, p2 - p1);
  grad[1] = inv6v * cross(p2 - p0, p3 - p0);
  grad[2] = inv6v * cross(p3 - p0, p1 - p0);
  grad[3] = inv6v * cross(p1 - p0, p2 - p0);
}

void tet_B_matrix(const Vec3 grad[4], double B[6][12]) {
  for (int i = 0; i < 6; ++i)
    for (int j = 0; j < 12; ++j) B[i][j] = 0.0;
  for (int a = 0; a < 4; ++a) {
    const double gx = grad[a].x, gy = grad[a].y, gz = grad[a].z;
    const int c0 = 3 * a, c1 = 3 * a + 1, c2 = 3 * a + 2;
    B[0][c0] = gx;
    B[1][c1] = gy;
    B[2][c2] = gz;
    B[3][c1] = gz;  B[3][c2] = gy;  // yz
    B[4][c0] = gz;  B[4][c2] = gx;  // xz
    B[5][c0] = gy;  B[5][c1] = gx;  // xy
  }
}

void iso_D_matrix(double E, double nu, double D[6][6]) {
  for (int i = 0; i < 6; ++i)
    for (int j = 0; j < 6; ++j) D[i][j] = 0.0;
  const double lam = E * nu / ((1.0 + nu) * (1.0 - 2.0 * nu));
  const double mu = E / (2.0 * (1.0 + nu));
  for (int i = 0; i < 3; ++i)
    for (int j = 0; j < 3; ++j) D[i][j] = lam + (i == j ? 2.0 * mu : 0.0);
  for (int i = 3; i < 6; ++i) D[i][i] = mu;
}

namespace {

// D * v for a 6-vector.
void dmul(const double D[6][6], const double v[6], double out[6]) {
  for (int i = 0; i < 6; ++i) {
    double s = 0.0;
    for (int j = 0; j < 6; ++j) s += D[i][j] * v[j];
    out[i] = s;
  }
}

}  // namespace

FemProblem fem_assemble(const Mesh& m, const std::vector<double>& E_cell,
                        const std::vector<double>& nu_cell,
                        const std::vector<double>& eps0_cell) {
  const int nn = static_cast<int>(m.nodes.size());
  const int ndof = 3 * nn;

  MeshTopology topo;
  topo.build(m);

  // Build the scalar CSR sparsity pattern: row block for node a covers
  // columns for {a} U node_adj[a] (every tet edge connects all 4 vertex
  // pairs, so plain edge adjacency already captures full element coupling).
  FemProblem prob;
  prob.K.n = ndof;
  prob.K.ptr.assign(ndof + 1, 0);
  std::vector<std::vector<int>> neighbors(nn);
  for (int a = 0; a < nn; ++a) {
    std::vector<int> nb = topo.node_adj[a];
    nb.push_back(a);
    std::sort(nb.begin(), nb.end());
    nb.erase(std::unique(nb.begin(), nb.end()), nb.end());
    neighbors[a] = std::move(nb);
  }
  for (int a = 0; a < nn; ++a) {
    const int rowblock = static_cast<int>(neighbors[a].size()) * 3;
    for (int comp = 0; comp < 3; ++comp) prob.K.ptr[3 * a + comp + 1] = rowblock;
  }
  for (int i = 0; i < ndof; ++i) prob.K.ptr[i + 1] += prob.K.ptr[i];
  const int nnz = prob.K.ptr[ndof];
  prob.K.col.assign(nnz, 0);
  prob.K.val.assign(nnz, 0.0);
  for (int a = 0; a < nn; ++a) {
    for (int comp = 0; comp < 3; ++comp) {
      const int row = 3 * a + comp;
      int k = prob.K.ptr[row];
      for (int b : neighbors[a])
        for (int d = 0; d < 3; ++d) prob.K.col[k++] = 3 * b + d;
    }
  }

  prob.f.assign(ndof, 0.0);
  prob.dof_fixed.assign(ndof, 0);
  prob.dof_val.assign(ndof, 0.0);

  // Element assembly.
  const int nc = static_cast<int>(m.cells.size());
  for (int c = 0; c < nc; ++c) {
    const auto& cell = m.cells[c];
    Vec3 p[4];
    for (int a = 0; a < 4; ++a) p[a] = m.nodes[cell[a]];
    double vol;
    Vec3 grad[4];
    tet_shape_grads(p, vol, grad);
    // Consistency check on the closed-form gradients (see fem.hpp doc).
    const Vec3 gsum = grad[0] + grad[1] + grad[2] + grad[3];
    assert(std::fabs(gsum.x) < 1e-9 && std::fabs(gsum.y) < 1e-9 &&
          std::fabs(gsum.z) < 1e-9);

    double B[6][12];
    tet_B_matrix(grad, B);
    double D[6][6];
    iso_D_matrix(E_cell[c], nu_cell[c], D);

    // ke = vol * B^T D B  (12x12); accumulate D*B columnwise first.
    double DB[6][12];
    for (int i = 0; i < 6; ++i)
      for (int j = 0; j < 12; ++j) {
        double s = 0.0;
        for (int k = 0; k < 6; ++k) s += D[i][k] * B[k][j];
        DB[i][j] = s;
      }
    double ke[12][12];
    for (int i = 0; i < 12; ++i)
      for (int j = 0; j < 12; ++j) {
        double s = 0.0;
        for (int k = 0; k < 6; ++k) s += B[k][i] * DB[k][j];
        ke[i][j] = vol * s;
      }

    // fe = vol * B^T D eps0
    const double* eps0 = &eps0_cell[c * 6];
    double Deps0[6];
    dmul(D, eps0, Deps0);
    double fe[12];
    for (int i = 0; i < 12; ++i) {
      double s = 0.0;
      for (int k = 0; k < 6; ++k) s += B[k][i] * Deps0[k];
      fe[i] = vol * s;
    }

    // Scatter into global K/f.
    int gdof[12];
    for (int a = 0; a < 4; ++a)
      for (int comp = 0; comp < 3; ++comp) gdof[3 * a + comp] = 3 * cell[a] + comp;
    for (int i = 0; i < 12; ++i) {
      prob.f[gdof[i]] += fe[i];
      for (int j = 0; j < 12; ++j) {
        const int slot = prob.K.find(gdof[i], gdof[j]);
        prob.K.val[slot] += ke[i][j];
      }
    }
  }

  return prob;
}

void fem_apply_bc(FemProblem& p) {
  const int n = p.K.n;
  for (int i = 0; i < n; ++i) {
    if (!p.dof_fixed[i]) continue;
    const double val = p.dof_val[i];
    // Eliminate column i's contribution from every other row's rhs, then
    // zero the off-diagonal column entries (symmetric elimination).
    for (int row = 0; row < n; ++row) {
      if (row == i) continue;
      const int slot = p.K.find(row, i);
      if (slot < 0) continue;
      p.f[row] -= p.K.val[slot] * val;
      p.K.val[slot] = 0.0;
    }
    // Turn row i into the identity row.
    for (int k = p.K.ptr[i]; k < p.K.ptr[i + 1]; ++k)
      p.K.val[k] = (p.K.col[k] == i) ? 1.0 : 0.0;
    p.f[i] = val;
  }
}

SolveResult fem_solve(const FemProblem& p, std::vector<double>& u,
                      double rtol, int maxit) {
  return cg_ilu0(p.K, p.f, u, rtol, maxit);
}

std::vector<double> fem_element_stress(const Mesh& m,
                                       const std::vector<double>& u,
                                       const std::vector<double>& E_cell,
                                       const std::vector<double>& nu_cell,
                                       const std::vector<double>& eps0_cell) {
  const int nc = static_cast<int>(m.cells.size());
  std::vector<double> stress(nc * 6, 0.0);
  for (int c = 0; c < nc; ++c) {
    const auto& cell = m.cells[c];
    Vec3 p[4];
    for (int a = 0; a < 4; ++a) p[a] = m.nodes[cell[a]];
    double vol;
    Vec3 grad[4];
    tet_shape_grads(p, vol, grad);
    double B[6][12];
    tet_B_matrix(grad, B);
    double D[6][6];
    iso_D_matrix(E_cell[c], nu_cell[c], D);

    double ul[12];
    for (int a = 0; a < 4; ++a)
      for (int comp = 0; comp < 3; ++comp) ul[3 * a + comp] = u[3 * cell[a] + comp];

    double Bu[6];
    for (int i = 0; i < 6; ++i) {
      double s = 0.0;
      for (int j = 0; j < 12; ++j) s += B[i][j] * ul[j];
      Bu[i] = s;
    }
    const double* eps0 = &eps0_cell[c * 6];
    double strain_mech[6];
    for (int i = 0; i < 6; ++i) strain_mech[i] = Bu[i] - eps0[i];
    double sigma[6];
    dmul(D, strain_mech, sigma);
    for (int i = 0; i < 6; ++i) stress[c * 6 + i] = sigma[i];
  }
  return stress;
}

}  // namespace cp
