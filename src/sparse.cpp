#include "cprocess/sparse.hpp"

#include <algorithm>
#include <cmath>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace cp {

// ---------------------------------------------------------------------------
// Parallel BLAS-1 / SpMV primitives
//
// Each row of a CSR SpMV is independent, dot products reduce, and the vector
// updates (AXPY) are elementwise — all trivially data-parallel. OpenMP is a
// no-op when the library is built without it, so these stay correct serially.
// ---------------------------------------------------------------------------

void CSR::mul(const std::vector<double>& x, std::vector<double>& y) const {
  y.assign(n, 0.0);
#pragma omp parallel for schedule(static)
  for (int i = 0; i < n; ++i) {
    double s = 0;
    for (int k = ptr[i]; k < ptr[i + 1]; ++k) s += val[k] * x[col[k]];
    y[i] = s;
  }
}

int CSR::find(int row, int c) const {
  int lo = ptr[row], hi = ptr[row + 1];
  while (lo < hi) {
    const int mid = (lo + hi) / 2;
    if (col[mid] < c)
      lo = mid + 1;
    else
      hi = mid;
  }
  return (lo < ptr[row + 1] && col[lo] == c) ? lo : -1;
}

namespace {

double dotv(const std::vector<double>& a, const std::vector<double>& b) {
  const int n = static_cast<int>(a.size());
  double s = 0;
#pragma omp parallel for reduction(+ : s) schedule(static)
  for (int i = 0; i < n; ++i) s += a[i] * b[i];
  return s;
}

// y <- y + alpha*x
void axpy(double alpha, const std::vector<double>& x, std::vector<double>& y) {
  const int n = static_cast<int>(x.size());
#pragma omp parallel for schedule(static)
  for (int i = 0; i < n; ++i) y[i] += alpha * x[i];
}

std::vector<double> inv_diag(const CSR& A) {
  std::vector<double> d(A.n, 1.0);
#pragma omp parallel for schedule(static)
  for (int i = 0; i < A.n; ++i) {
    const int k = A.find(i, i);
    if (k >= 0 && A.val[k] != 0.0) d[i] = 1.0 / A.val[k];
  }
  return d;
}

}  // namespace

// ---------------------------------------------------------------------------
// ILU(0) with level-scheduled parallel triangular solves
// ---------------------------------------------------------------------------

void ILU0::factor(const CSR& A) {
  const int n = A.n;
  lu = A;  // copy pattern + values; we overwrite values in place
  diag.assign(n, -1);
  for (int i = 0; i < n; ++i) diag[i] = lu.find(i, i);

  // Standard IKJ incomplete LU with no fill-in (values only where A is nonzero).
  // This setup phase is inherently sequential (row i depends on rows < i).
  for (int i = 0; i < n; ++i) {
    for (int kk = lu.ptr[i]; kk < lu.ptr[i + 1]; ++kk) {
      const int k = lu.col[kk];
      if (k >= i) break;  // columns are sorted; L part is k < i
      const int dk = diag[k];
      if (dk < 0 || lu.val[dk] == 0.0) continue;
      const double lik = lu.val[kk] / lu.val[dk];
      lu.val[kk] = lik;
      // Update the rest of row i: a_ij -= lik * u_kj for j > k that exist in i.
      int jj = kk + 1;
      for (int pk = dk + 1; pk < lu.ptr[k + 1] && jj < lu.ptr[i + 1]; ++pk) {
        const int j = lu.col[pk];
        while (jj < lu.ptr[i + 1] && lu.col[jj] < j) ++jj;
        if (jj < lu.ptr[i + 1] && lu.col[jj] == j)
          lu.val[jj] -= lik * lu.val[pk];
      }
    }
  }

  // Level scheduling. lower solve: row i depends on rows j<i with L_ij!=0.
  // upper solve: row i depends on rows j>i with U_ij!=0. A row's level is one
  // past the max level of its dependencies; rows in one level are independent.
  std::vector<int> lev_lo(n, 0), lev_up(n, 0);
  int nlo = 0, nup = 0;
  for (int i = 0; i < n; ++i) {
    int lv = 0;
    for (int kk = lu.ptr[i]; kk < lu.ptr[i + 1]; ++kk) {
      const int j = lu.col[kk];
      if (j < i) lv = std::max(lv, lev_lo[j] + 1);
    }
    lev_lo[i] = lv;
    nlo = std::max(nlo, lv + 1);
  }
  for (int i = n - 1; i >= 0; --i) {
    int lv = 0;
    for (int kk = lu.ptr[i]; kk < lu.ptr[i + 1]; ++kk) {
      const int j = lu.col[kk];
      if (j > i) lv = std::max(lv, lev_up[j] + 1);
    }
    lev_up[i] = lv;
    nup = std::max(nup, lv + 1);
  }
  lvl_lo.assign(nlo, {});
  lvl_up.assign(nup, {});
  for (int i = 0; i < n; ++i) lvl_lo[lev_lo[i]].push_back(i);
  for (int i = 0; i < n; ++i) lvl_up[lev_up[i]].push_back(i);
}

void ILU0::apply(const std::vector<double>& x, std::vector<double>& y) const {
  const int n = lu.n;
  y.resize(n);
  // Forward solve L z = x (L is unit lower; stored strictly-lower part).
  for (const auto& level : lvl_lo) {
    const int m = static_cast<int>(level.size());
#pragma omp parallel for schedule(static)
    for (int t = 0; t < m; ++t) {
      const int i = level[t];
      double s = x[i];
      for (int kk = lu.ptr[i]; kk < lu.ptr[i + 1]; ++kk) {
        const int j = lu.col[kk];
        if (j < i) s -= lu.val[kk] * y[j];
      }
      y[i] = s;  // unit diagonal
    }
  }
  // Backward solve U y = z (U includes the diagonal).
  for (const auto& level : lvl_up) {
    const int m = static_cast<int>(level.size());
#pragma omp parallel for schedule(static)
    for (int t = 0; t < m; ++t) {
      const int i = level[t];
      double s = y[i];
      double dii = 1.0;
      for (int kk = lu.ptr[i]; kk < lu.ptr[i + 1]; ++kk) {
        const int j = lu.col[kk];
        if (j > i) s -= lu.val[kk] * y[j];
        else if (j == i) dii = lu.val[kk];
      }
      y[i] = (dii != 0.0) ? s / dii : s;
    }
  }
}

// ---------------------------------------------------------------------------
// Preconditioned Krylov cores
// ---------------------------------------------------------------------------

SolveResult cg(const CSR& A, const std::vector<double>& b,
               std::vector<double>& x, double rtol, int maxit,
               const Precond& psolve) {
  SolveResult res;
  const int n = A.n;
  x.resize(n, 0.0);
  const double bnorm = std::sqrt(dotv(b, b));
  if (bnorm == 0.0) {
    x.assign(n, 0.0);
    res.converged = true;
    return res;
  }
  std::vector<double> r(n), z(n), p(n), q(n);
  A.mul(x, q);
#pragma omp parallel for schedule(static)
  for (int i = 0; i < n; ++i) r[i] = b[i] - q[i];
  psolve(r, z);
  p = z;
  double rz = dotv(r, z);
  for (int it = 0; it < maxit; ++it) {
    A.mul(p, q);
    const double pq = dotv(p, q);
    if (pq == 0.0) break;
    const double alpha = rz / pq;
    axpy(alpha, p, x);
    axpy(-alpha, q, r);
    const double rn = std::sqrt(dotv(r, r));
    res.iters = it + 1;
    res.resid = rn / bnorm;
    if (rn <= rtol * bnorm) {
      res.converged = true;
      return res;
    }
    psolve(r, z);
    const double rz2 = dotv(r, z);
    const double beta = rz2 / rz;
    rz = rz2;
#pragma omp parallel for schedule(static)
    for (int i = 0; i < n; ++i) p[i] = z[i] + beta * p[i];
  }
  return res;
}

SolveResult bicgstab(const CSR& A, const std::vector<double>& b,
                     std::vector<double>& x, double rtol, int maxit,
                     const Precond& psolve) {
  SolveResult res;
  const int n = A.n;
  x.resize(n, 0.0);
  const double bnorm = std::sqrt(dotv(b, b));
  if (bnorm == 0.0) {
    x.assign(n, 0.0);
    res.converged = true;
    return res;
  }
  std::vector<double> r(n), r0(n), p(n, 0.0), v(n, 0.0), s(n), t(n), ph(n), sh(n);
  A.mul(x, t);
#pragma omp parallel for schedule(static)
  for (int i = 0; i < n; ++i) r[i] = b[i] - t[i];
  r0 = r;
  double rho = 1, alpha = 1, omega = 1;
  for (int it = 0; it < maxit; ++it) {
    const double rho2 = dotv(r0, r);
    if (rho2 == 0.0) break;
    if (it == 0) {
      p = r;
    } else {
      const double beta = (rho2 / rho) * (alpha / omega);
#pragma omp parallel for schedule(static)
      for (int i = 0; i < n; ++i) p[i] = r[i] + beta * (p[i] - omega * v[i]);
    }
    rho = rho2;
    psolve(p, ph);
    A.mul(ph, v);
    const double r0v = dotv(r0, v);
    if (r0v == 0.0) break;
    alpha = rho / r0v;
#pragma omp parallel for schedule(static)
    for (int i = 0; i < n; ++i) s[i] = r[i] - alpha * v[i];
    double sn = std::sqrt(dotv(s, s));
    res.iters = it + 1;
    if (sn <= rtol * bnorm) {
      axpy(alpha, ph, x);
      res.resid = sn / bnorm;
      res.converged = true;
      return res;
    }
    psolve(s, sh);
    A.mul(sh, t);
    const double tt = dotv(t, t);
    if (tt == 0.0) break;
    omega = dotv(t, s) / tt;
#pragma omp parallel for schedule(static)
    for (int i = 0; i < n; ++i) x[i] += alpha * ph[i] + omega * sh[i];
#pragma omp parallel for schedule(static)
    for (int i = 0; i < n; ++i) r[i] = s[i] - omega * t[i];
    const double rn = std::sqrt(dotv(r, r));
    res.resid = rn / bnorm;
    if (rn <= rtol * bnorm) {
      res.converged = true;
      return res;
    }
    if (omega == 0.0) break;
  }
  return res;
}

// ---------------------------------------------------------------------------
// Public convenience wrappers
// ---------------------------------------------------------------------------

SolveResult cg_jacobi(const CSR& A, const std::vector<double>& b,
                      std::vector<double>& x, double rtol, int maxit) {
  const std::vector<double> M = inv_diag(A);
  Precond psolve = [&](const std::vector<double>& in, std::vector<double>& out) {
    const int n = A.n;
    out.resize(n);
#pragma omp parallel for schedule(static)
    for (int i = 0; i < n; ++i) out[i] = M[i] * in[i];
  };
  return cg(A, b, x, rtol, maxit, psolve);
}

SolveResult bicgstab_jacobi(const CSR& A, const std::vector<double>& b,
                            std::vector<double>& x, double rtol, int maxit) {
  const std::vector<double> M = inv_diag(A);
  Precond psolve = [&](const std::vector<double>& in, std::vector<double>& out) {
    const int n = A.n;
    out.resize(n);
#pragma omp parallel for schedule(static)
    for (int i = 0; i < n; ++i) out[i] = M[i] * in[i];
  };
  return bicgstab(A, b, x, rtol, maxit, psolve);
}

SolveResult cg_ilu0(const CSR& A, const std::vector<double>& b,
                    std::vector<double>& x, double rtol, int maxit) {
  ILU0 ilu;
  ilu.factor(A);
  Precond psolve = [&](const std::vector<double>& in, std::vector<double>& out) {
    ilu.apply(in, out);
  };
  return cg(A, b, x, rtol, maxit, psolve);
}

SolveResult bicgstab_ilu0(const CSR& A, const std::vector<double>& b,
                          std::vector<double>& x, double rtol, int maxit) {
  ILU0 ilu;
  ilu.factor(A);
  Precond psolve = [&](const std::vector<double>& in, std::vector<double>& out) {
    ilu.apply(in, out);
  };
  return bicgstab(A, b, x, rtol, maxit, psolve);
}

// ---------------------------------------------------------------------------
// GMRES(restart) with right Jacobi preconditioning
//
// Algorithm: Saad & Schultz 1986 with Givens-rotation QR update.
// Right preconditioning (solve A*M^{-1}*(M*x) = b, then x <- M^{-1}*y):
//   - true residual is ||b - A*x|| at every restart
//   - preconditioner affects the search direction, not the residual norm
//
// Storage:
//   V[j][i] : j-th Arnoldi basis vector, component i  [(restart+1) x n]
//   H[i*rs+j]: Hessenberg element (i,j), row-major    [(restart+1)*restart]
//   c[j],s[j]: Givens rotation applied at step j
//   g[j]     : RHS of the projected least-squares problem (rotated)
// ---------------------------------------------------------------------------
SolveResult gmres_jacobi(const CSR& A, const std::vector<double>& b,
                         std::vector<double>& x, double rtol, int maxit,
                         int restart) {
  SolveResult res;
  const int n = A.n;
  x.resize(n, 0.0);
  const double bnorm = std::sqrt(dotv(b, b));
  if (bnorm == 0.0) {
    x.assign(n, 0.0);
    res.converged = true;
    return res;
  }

  const std::vector<double> M = inv_diag(A);
  std::vector<double> r(n), t(n);

  for (int total = 0; total < maxit; ) {
    // r = b - A*x (true residual at restart boundary)
    A.mul(x, t);
#pragma omp parallel for schedule(static)
    for (int i = 0; i < n; ++i) r[i] = b[i] - t[i];
    const double beta = std::sqrt(dotv(r, r));
    res.resid = beta / bnorm;
    if (res.resid <= rtol) { res.converged = true; return res; }
    if (beta < 1e-14)       { res.converged = true; return res; }

    // Arnoldi basis V[0..restart], Hessenberg H[(restart+1)*restart]
    const int rs = restart;
    std::vector<std::vector<double>> V(rs + 1, std::vector<double>(n, 0.0));
    std::vector<double> H((rs + 1) * rs, 0.0);
    std::vector<double> g(rs + 1, 0.0);
    std::vector<double> c(rs, 0.0), s(rs, 0.0);

#pragma omp parallel for schedule(static)
    for (int i = 0; i < n; ++i) V[0][i] = r[i] / beta;
    g[0] = beta;

    int m = 0;  // columns of H built so far
    double hnext = 0.0;
    for (; m < rs && total < maxit; ++m, ++total) {
      // w = A * M^{-1} * V[m]  (right-preconditioned matrix-vector product)
#pragma omp parallel for schedule(static)
      for (int i = 0; i < n; ++i) t[i] = M[i] * V[m][i];
      A.mul(t, r);  // r = w

      // Modified Gram-Schmidt: orthogonalize w against V[0..m]
      for (int j = 0; j <= m; ++j) {
        const double h = dotv(r, V[j]);
        H[j * rs + m] = h;
        axpy(-h, V[j], r);
      }
      hnext = std::sqrt(dotv(r, r));
      H[(m + 1) * rs + m] = hnext;

      if (hnext > 1e-12)
#pragma omp parallel for schedule(static)
        for (int i = 0; i < n; ++i) V[m + 1][i] = r[i] / hnext;

      // Apply previous Givens rotations to column m of H
      for (int j = 0; j < m; ++j) {
        const double h0 = H[j * rs + m], h1 = H[(j + 1) * rs + m];
        H[j * rs + m]       =  c[j] * h0 + s[j] * h1;
        H[(j + 1) * rs + m] = -s[j] * h0 + c[j] * h1;
      }

      // Compute new Givens rotation to zero H[(m+1),m]
      const double a = H[m * rs + m], bh = H[(m + 1) * rs + m];
      const double rho = std::sqrt(a * a + bh * bh);
      if (rho < 1e-14) { ++m; break; }
      c[m] = a / rho;  s[m] = bh / rho;
      H[m * rs + m]       = rho;
      H[(m + 1) * rs + m] = 0.0;

      // Update g: apply the new rotation
      g[m + 1] = -s[m] * g[m];
      g[m]     =  c[m] * g[m];

      res.resid = std::fabs(g[m + 1]) / bnorm;
      res.iters = total + 1;
      if (res.resid <= rtol) { ++m; break; }
    }

    // Back-solve upper triangular R (m x m) stored in H
    std::vector<double> y(m, 0.0);
    for (int i = m - 1; i >= 0; --i) {
      y[i] = g[i];
      for (int j = i + 1; j < m; ++j) y[i] -= H[i * rs + j] * y[j];
      y[i] /= H[i * rs + i];
    }

    // Update x: x += M^{-1} * (V_m * y)  (undo right preconditioning)
#pragma omp parallel for schedule(static)
    for (int i = 0; i < n; ++i) {
      double acc = 0.0;
      for (int j = 0; j < m; ++j) acc += V[j][i] * y[j];
      x[i] += M[i] * acc;
    }

    if (res.resid <= rtol) { res.converged = true; return res; }
    if (hnext < 1e-12) break;  // exact solution found mid-restart
  }

  // Recompute true residual after final update
  A.mul(x, t);
#pragma omp parallel for schedule(static)
  for (int i = 0; i < n; ++i) r[i] = b[i] - t[i];
  res.resid = std::sqrt(dotv(r, r)) / bnorm;
  res.converged = (res.resid <= rtol);
  return res;
}

}  // namespace cp
