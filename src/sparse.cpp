#include "cprocess/sparse.hpp"

#include <algorithm>
#include <cmath>

namespace cp {

void CSR::mul(const std::vector<double>& x, std::vector<double>& y) const {
  y.assign(n, 0.0);
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

static double dotv(const std::vector<double>& a, const std::vector<double>& b) {
  double s = 0;
  for (std::size_t i = 0; i < a.size(); ++i) s += a[i] * b[i];
  return s;
}

static std::vector<double> inv_diag(const CSR& A) {
  std::vector<double> d(A.n, 1.0);
  for (int i = 0; i < A.n; ++i) {
    const int k = A.find(i, i);
    if (k >= 0 && A.val[k] != 0.0) d[i] = 1.0 / A.val[k];
  }
  return d;
}

SolveResult cg_jacobi(const CSR& A, const std::vector<double>& b,
                      std::vector<double>& x, double rtol, int maxit) {
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
  std::vector<double> r(n), z(n), p(n), q(n);
  A.mul(x, q);
  for (int i = 0; i < n; ++i) r[i] = b[i] - q[i];
  for (int i = 0; i < n; ++i) z[i] = M[i] * r[i];
  p = z;
  double rz = dotv(r, z);
  for (int it = 0; it < maxit; ++it) {
    A.mul(p, q);
    const double pq = dotv(p, q);
    if (pq == 0.0) break;
    const double alpha = rz / pq;
    for (int i = 0; i < n; ++i) x[i] += alpha * p[i];
    for (int i = 0; i < n; ++i) r[i] -= alpha * q[i];
    const double rn = std::sqrt(dotv(r, r));
    res.iters = it + 1;
    res.resid = rn / bnorm;
    if (rn <= rtol * bnorm) {
      res.converged = true;
      return res;
    }
    for (int i = 0; i < n; ++i) z[i] = M[i] * r[i];
    const double rz2 = dotv(r, z);
    const double beta = rz2 / rz;
    rz = rz2;
    for (int i = 0; i < n; ++i) p[i] = z[i] + beta * p[i];
  }
  return res;
}

SolveResult bicgstab_jacobi(const CSR& A, const std::vector<double>& b,
                            std::vector<double>& x, double rtol, int maxit) {
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
  std::vector<double> r(n), r0(n), p(n, 0.0), v(n, 0.0), s(n), t(n), ph(n), sh(n);
  A.mul(x, t);
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
      for (int i = 0; i < n; ++i) p[i] = r[i] + beta * (p[i] - omega * v[i]);
    }
    rho = rho2;
    for (int i = 0; i < n; ++i) ph[i] = M[i] * p[i];
    A.mul(ph, v);
    const double r0v = dotv(r0, v);
    if (r0v == 0.0) break;
    alpha = rho / r0v;
    for (int i = 0; i < n; ++i) s[i] = r[i] - alpha * v[i];
    double sn = std::sqrt(dotv(s, s));
    res.iters = it + 1;
    if (sn <= rtol * bnorm) {
      for (int i = 0; i < n; ++i) x[i] += alpha * ph[i];
      res.resid = sn / bnorm;
      res.converged = true;
      return res;
    }
    for (int i = 0; i < n; ++i) sh[i] = M[i] * s[i];
    A.mul(sh, t);
    const double tt = dotv(t, t);
    if (tt == 0.0) break;
    omega = dotv(t, s) / tt;
    for (int i = 0; i < n; ++i) x[i] += alpha * ph[i] + omega * sh[i];
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

}  // namespace cp
