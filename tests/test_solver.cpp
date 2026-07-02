#include <cmath>
#include <cstdio>
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif

#include "cprocess/sparse.hpp"
#include "test_util.hpp"

using namespace cp;

// 1D Laplacian (tridiagonal SPD), manufactured solution.
static CSR laplacian1d(int n) {
  CSR a;
  a.n = n;
  a.ptr.push_back(0);
  for (int i = 0; i < n; ++i) {
    if (i > 0) { a.col.push_back(i - 1); a.val.push_back(-1.0); }
    a.col.push_back(i); a.val.push_back(2.0);
    if (i < n - 1) { a.col.push_back(i + 1); a.val.push_back(-1.0); }
    a.ptr.push_back(static_cast<int>(a.col.size()));
  }
  return a;
}

// 2D 5-point Laplacian on an nx*ny grid (SPD), row-major node ordering.
// Wide bandwidth gives level scheduling real parallelism (unlike the
// tridiagonal 1D case, whose rows form a strictly sequential chain).
static CSR laplacian2d(int nx, int ny) {
  CSR a;
  a.n = nx * ny;
  a.ptr.push_back(0);
  for (int iy = 0; iy < ny; ++iy) {
    for (int ix = 0; ix < nx; ++ix) {
      const int i = iy * nx + ix;
      if (iy > 0) { a.col.push_back(i - nx); a.val.push_back(-1.0); }
      if (ix > 0) { a.col.push_back(i - 1); a.val.push_back(-1.0); }
      a.col.push_back(i); a.val.push_back(4.0);
      if (ix < nx - 1) { a.col.push_back(i + 1); a.val.push_back(-1.0); }
      if (iy < ny - 1) { a.col.push_back(i + nx); a.val.push_back(-1.0); }
      a.ptr.push_back(static_cast<int>(a.col.size()));
    }
  }
  return a;
}

int main() {
  const int n = 200;
  CSR a = laplacian1d(n);

  std::vector<double> xtrue(n), b(n), x;
  for (int i = 0; i < n; ++i) xtrue[i] = std::sin(0.1 * i) + 0.5;
  a.mul(xtrue, b);

  SolveResult r = cg_jacobi(a, b, x, 1e-12, 2000);
  CHECK(r.converged);
  for (int i = 0; i < n; ++i) CHECK_NEAR(x[i], xtrue[i], 1e-7);
  std::printf("cg: %d iters, resid %.2e\n", r.iters, r.resid);

  x.assign(n, 0.0);
  r = bicgstab_jacobi(a, b, x, 1e-12, 2000);
  CHECK(r.converged);
  for (int i = 0; i < n; ++i) CHECK_NEAR(x[i], xtrue[i], 1e-6);
  std::printf("bicgstab: %d iters, resid %.2e\n", r.iters, r.resid);

  // GMRES: use n=50 so GMRES(30) converges within the restart budget.
  // (1D Laplacian has κ ~ n^2; for n=200 GMRES(30) needs too many restarts.)
  {
    const int ng = 50;
    CSR ag = laplacian1d(ng);
    std::vector<double> xg, bg(ng), xtg(ng);
    for (int i = 0; i < ng; ++i) xtg[i] = std::sin(0.1*i) + 0.5;
    ag.mul(xtg, bg);
    SolveResult rg = gmres_jacobi(ag, bg, xg, 1e-12, 2000, 30);
    CHECK(rg.converged);
    for (int i = 0; i < ng; ++i) CHECK_NEAR(xg[i], xtg[i], 1e-6);
    std::printf("gmres(30): %d iters, resid %.2e\n", rg.iters, rg.resid);
  }

  // ILU(0)-preconditioned CG: same manufactured solution, should converge in
  // fewer iterations than plain Jacobi (ILU(0) is exact for a tridiagonal M).
  x.assign(n, 0.0);
  SolveResult ri = cg_ilu0(a, b, x, 1e-12, 2000);
  CHECK(ri.converged);
  for (int i = 0; i < n; ++i) CHECK_NEAR(x[i], xtrue[i], 1e-7);
  std::printf("cg_ilu0: %d iters, resid %.2e\n", ri.iters, ri.resid);
  CHECK(ri.iters <= r.iters);  // ILU(0) never worse than Jacobi here

  x.assign(n, 0.0);
  ri = bicgstab_ilu0(a, b, x, 1e-12, 2000);
  CHECK(ri.converged);
  for (int i = 0; i < n; ++i) CHECK_NEAR(x[i], xtrue[i], 1e-6);
  std::printf("bicgstab_ilu0: %d iters, resid %.2e\n", ri.iters, ri.resid);

  // ILU(0) factor + apply directly solves the tridiagonal system exactly
  // (no fill-in is dropped for a tridiagonal matrix), so M^{-1} b == A^{-1} b.
  {
    ILU0 ilu;
    ilu.factor(a);
    std::vector<double> y;
    ilu.apply(b, y);
    for (int i = 0; i < n; ++i) CHECK_NEAR(y[i], xtrue[i], 1e-9);
    std::printf("ilu0 exact for tridiagonal: ok\n");
  }

  // Zero RHS edge case.
  std::vector<double> zb(n, 0.0);
  r = cg_jacobi(a, zb, x, 1e-12, 100);
  CHECK(r.converged);
  for (int i = 0; i < n; ++i) CHECK(x[i] == 0.0);

  // Level-parallel ILU(0) factorization: bit-identical results at thread
  // counts 1 and 4 on a 2D 5-point Laplacian, whose bandwidth gives level
  // scheduling real parallelism.
  {
    const int nx = 50, ny = 50;
    CSR a2 = laplacian2d(nx, ny);

#ifdef _OPENMP
    omp_set_num_threads(1);
#endif
    ILU0 ilu1;
    ilu1.factor(a2);

#ifdef _OPENMP
    omp_set_num_threads(4);
#endif
    ILU0 ilu4;
    ilu4.factor(a2);

    CHECK(ilu1.lu.val.size() == ilu4.lu.val.size());
    for (size_t k = 0; k < ilu1.lu.val.size(); ++k)
      CHECK(ilu1.lu.val[k] == ilu4.lu.val[k]);  // bit-identical
    std::printf("ilu0 factor: bit-identical across thread counts\n");

    // Convergence test on the same 2D Laplacian with a manufactured solution.
    const int n2 = nx * ny;
    std::vector<double> xtrue2(n2), b2, x2;
    for (int i = 0; i < n2; ++i) xtrue2[i] = std::sin(0.1 * i) + 0.5;
    a2.mul(xtrue2, b2);
    SolveResult r2 = cg_ilu0(a2, b2, x2, 1e-10, 500);
    CHECK(r2.converged);
    double maxerr = 0.0;
    for (int i = 0; i < n2; ++i)
      maxerr = std::max(maxerr, std::fabs(x2[i] - xtrue2[i]));
    CHECK(maxerr < 1e-6);
    std::printf("cg_ilu0 (2D laplacian, 4 threads): %d iters, resid %.2e, maxerr %.2e\n",
                r2.iters, r2.resid, maxerr);

#ifdef _OPENMP
    omp_set_num_threads(1);
#endif
  }

  std::printf("solver tests passed\n");
  return 0;
}
