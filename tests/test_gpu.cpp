// Only compiled/registered when CPROCESS_GPU is enabled (see CMakeLists.txt).
// Verifies the OpenMP target-offload Jacobi-CG (cg_jacobi_gpu) against the
// existing host cg_jacobi on a 2D 5-point Laplacian. On a development box
// without an offload device, the `target` regions execute on the host per
// the OpenMP spec, so this test validates correctness even there.
#include <cmath>
#include <cstdio>
#include <vector>

#include "cprocess/sparse.hpp"
#include "test_util.hpp"

using namespace cp;

// 2D 5-point Laplacian on an nx*ny grid (SPD), row-major node ordering.
// (Duplicated from tests/test_solver.cpp per the PA-4 spec.)
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
  const int nx = 50, ny = 50;
  CSR a = laplacian2d(nx, ny);
  const int n = nx * ny;

  std::vector<double> xtrue(n), b;
  for (int i = 0; i < n; ++i) xtrue[i] = std::sin(0.1 * i) + 0.5;
  a.mul(xtrue, b);

  // Correctness: GPU and CPU solutions agree to < 1e-9.
  std::vector<double> x_cpu, x_gpu;
  SolveResult r_cpu = cg_jacobi(a, b, x_cpu, 1e-12, 2000);
  SolveResult r_gpu = cg_jacobi_gpu(a, b, x_gpu, 1e-12, 2000);
  CHECK(r_cpu.converged);
  CHECK(r_gpu.converged);

  double maxdiff = 0.0;
  for (int i = 0; i < n; ++i)
    maxdiff = std::max(maxdiff, std::fabs(x_cpu[i] - x_gpu[i]));
  CHECK(maxdiff < 1e-9);
  std::printf("cg_jacobi_gpu vs cg_jacobi: maxdiff %.2e\n", maxdiff);

  // Convergence: iteration count within +/-20% of the CPU version.
  CHECK(r_gpu.iters <= static_cast<int>(1.2 * r_cpu.iters) + 1);
  CHECK(r_gpu.iters >= static_cast<int>(0.8 * r_cpu.iters) - 1);
  std::printf("cg_jacobi_gpu: %d iters (cpu %d), resid %.2e\n", r_gpu.iters,
              r_cpu.iters, r_gpu.resid);

  // Zero RHS edge case: x=0, converged=true.
  std::vector<double> zb(n, 0.0), xz;
  SolveResult rz = cg_jacobi_gpu(a, zb, xz, 1e-12, 100);
  CHECK(rz.converged);
  for (int i = 0; i < n; ++i) CHECK(xz[i] == 0.0);

  std::printf("gpu tests passed\n");
  return 0;
}
