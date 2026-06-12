#include <cmath>
#include <cstdio>
#include <vector>

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

  // Zero RHS edge case.
  std::vector<double> zb(n, 0.0);
  r = cg_jacobi(a, zb, x, 1e-12, 100);
  CHECK(r.converged);
  for (int i = 0; i < n; ++i) CHECK(x[i] == 0.0);

  std::printf("solver tests passed\n");
  return 0;
}
