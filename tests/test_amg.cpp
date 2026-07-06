#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif

#include "cprocess/sparse.hpp"
#include "test_util.hpp"

using namespace cp;

// n^3 grid, 7-point stencil, diagonal 6, neighbors -1, Dirichlet boundary
// (i.e. missing neighbors at the boundary are simply omitted — standard
// truncated-stencil Dirichlet Laplacian).
static CSR laplacian3d(int n) {
  CSR a;
  a.n = n * n * n;
  a.ptr.push_back(0);
  auto idx = [n](int x, int y, int z) { return (z * n + y) * n + x; };
  for (int z = 0; z < n; ++z) {
    for (int y = 0; y < n; ++y) {
      for (int x = 0; x < n; ++x) {
        std::vector<std::pair<int, double>> row;
        if (z > 0) row.emplace_back(idx(x, y, z - 1), -1.0);
        if (y > 0) row.emplace_back(idx(x, y - 1, z), -1.0);
        if (x > 0) row.emplace_back(idx(x - 1, y, z), -1.0);
        row.emplace_back(idx(x, y, z), 6.0);
        if (x < n - 1) row.emplace_back(idx(x + 1, y, z), -1.0);
        if (y < n - 1) row.emplace_back(idx(x, y + 1, z), -1.0);
        if (z < n - 1) row.emplace_back(idx(x, y, z + 1), -1.0);
        std::sort(row.begin(), row.end());
        for (auto& pr : row) { a.col.push_back(pr.first); a.val.push_back(pr.second); }
        a.ptr.push_back(static_cast<int>(a.col.size()));
      }
    }
  }
  return a;
}

int main() {
  // --- Test 1 & 2: iteration count and solution agreement, n=32 ---
  {
    const int n = 32;
    CSR A = laplacian3d(n);
    std::vector<double> b(A.n, 1.0), x_amg, x_ilu;

    auto t0 = std::chrono::steady_clock::now();
    SolveResult r_amg = cg_amg(A, b, x_amg, 1e-10, 500);
    auto t1 = std::chrono::steady_clock::now();
    SolveResult r_ilu = cg_ilu0(A, b, x_ilu, 1e-10, 500);
    auto t2 = std::chrono::steady_clock::now();

    std::printf("AMG: setup+solve time = %.4f s, iters = %d\n",
                std::chrono::duration<double>(t1 - t0).count(), r_amg.iters);
    std::printf("ILU0: solve time = %.4f s, iters = %d\n",
                std::chrono::duration<double>(t2 - t1).count(), r_ilu.iters);

    CHECK(r_amg.converged);
    CHECK(r_ilu.converged);
    CHECK(r_amg.iters < 40);
    CHECK(r_amg.iters * 2 < r_ilu.iters);

    double maxdiff = 0.0;
    for (int i = 0; i < A.n; ++i)
      maxdiff = std::max(maxdiff, std::fabs(x_amg[i] - x_ilu[i]));
    std::printf("max|x_amg - x_ilu| = %.3e\n", maxdiff);
    CHECK(maxdiff < 1e-8);

    // --- Test 4: aggregate count ---
    AMG amg;
    amg.setup(A);
    const int nagg = amg.n_aggregates();
    std::printf("n_aggregates = %d (n=%d)\n", nagg, A.n);
    CHECK(nagg >= A.n / 20);
    CHECK(nagg <= A.n / 3);
  }

  // --- Test 5: thread-count determinism, n=16 ---
  {
    const int n = 16;
    CSR A = laplacian3d(n);
    std::vector<double> b(A.n, 1.0);

#ifdef _OPENMP
    omp_set_num_threads(1);
#endif
    std::vector<double> x1;
    SolveResult r1 = cg_amg(A, b, x1, 1e-10, 500);
    CHECK(r1.converged);

#ifdef _OPENMP
    omp_set_num_threads(4);
#endif
    std::vector<double> x4;
    SolveResult r4 = cg_amg(A, b, x4, 1e-10, 500);
    CHECK(r4.converged);

    double num = 0.0, den = 0.0;
    for (int i = 0; i < A.n; ++i) {
      num = std::max(num, std::fabs(x1[i] - x4[i]));
      den = std::max(den, std::fabs(x1[i]));
    }
    const double relerr = (den > 0.0) ? num / den : num;
    std::printf("thread determinism relerr = %.3e\n", relerr);
    CHECK(relerr < 1e-9);
  }

  std::printf("test_amg: all checks passed\n");
  return 0;
}
