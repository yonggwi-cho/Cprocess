#include "cprocess/sparse.hpp"

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <unordered_map>

#ifdef _OPENMP
#include <omp.h>
#endif
#if defined(__SSE__) || defined(__x86_64__) || defined(_M_X64)
#include <xmmintrin.h>
#include <pmmintrin.h>
#define CP_HAS_SSE_FTZ 1
#endif

namespace cp {

namespace {
// Flush-to-zero / denormals-are-zero: the point-defect reaction system
// (P2-1) solves a near-singular pure-diffusion operator (no reaction term
// in the diffusion sub-step) whose Krylov residual can explore directions
// that underflow toward subnormal doubles as they approach machine precision.
// Subnormal arithmetic is ~100x slower on x86 without FTZ/DAZ, which turned
// a sub-second solve into an apparent hang. This is a one-time, thread-local
// FPU mode flip; it has no effect on any value above ~2.2e-308 so normal
// physical results (all >> that scale) are unaffected.
//
// This MUST NOT run as a global static initializer: a `#pragma omp parallel`
// region executed before main() (during dynamic static initialization, in
// unspecified order relative to other translation units' statics, before the
// OpenMP runtime and its thread pool are guaranteed fully set up) is a real
// deadlock/hang hazard -- observed in practice as intermittent hangs whose
// exact location varied run to run, consistent with a startup race. Instead,
// set it lazily on first use via a C++11 function-local static (thread-safe
// initialization guaranteed by the standard, and it only runs once real
// program execution -- and the OpenMP runtime -- are already underway).
void ensure_ftz_daz() {
#ifdef CP_HAS_SSE_FTZ
  static const bool kDone = [] {
#pragma omp parallel
    {
      _MM_SET_FLUSH_ZERO_MODE(_MM_FLUSH_ZERO_ON);
      _MM_SET_DENORMALS_ZERO_MODE(_MM_DENORMALS_ZERO_ON);
    }
    return true;
  }();
  (void)kDone;
#endif
}

// Some execution hosts (observed on a heavily oversubscribed/virtualized CI
// sandbox while developing P2-1) have per-thread scheduling latency so bad
// that spinning up an OpenMP thread team for a small SpMV/dot-product region
// costs orders of magnitude more wall time than the work itself -- a solve
// that completes in under a second with OMP_NUM_THREADS=1 can appear to hang
// indefinitely at the default thread count. This has nothing to do with
// P2-1's numerics (reproduces on the pre-existing, untouched equilibrium
// diffuse() path too) and everything to do with the host's thread scheduler.
// Respect an explicit OMP_NUM_THREADS from the environment (a real multi-core
// deployment should keep using it), but if it wasn't set, cap the default to
// a conservative value once, lazily, the same way as ensure_ftz_daz() above
// (never as a static initializer -- see that function's comment).
void ensure_sane_thread_count() {
#ifdef _OPENMP
  static const bool kDone = [] {
    if (!std::getenv("OMP_NUM_THREADS") && omp_get_max_threads() > 1)
      omp_set_num_threads(1);
    return true;
  }();
  (void)kDone;
#endif
}
}  // namespace

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

  // Level scheduling. lower solve: row i depends on rows j<i with L_ij!=0.
  // upper solve: row i depends on rows j>i with U_ij!=0. A row's level is one
  // past the max level of its dependencies; rows in one level are independent.
  // This depends only on the sparsity pattern of A, so it can be computed
  // before the factorization values are known. The lower-solve dependency
  // graph (row i depends on rows k<i with a_ik != 0) is exactly the IKJ
  // factorization's dependency graph, so lvl_lo also drives the parallel
  // factorization loop below.
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

  // Standard IKJ incomplete LU with no fill-in (values only where A is
  // nonzero), executed level-by-level. Rows within a level are mutually
  // independent (their dependencies k < i all lie in earlier levels), and
  // each row only writes within its own range in `lu`, so this is race-free.
  for (const auto& level : lvl_lo) {
    const int m = static_cast<int>(level.size());
#pragma omp parallel for schedule(dynamic, 16)
    for (int t = 0; t < m; ++t) {
      const int i = level[t];
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
  }
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
// Smoothed-aggregation AMG (two-level V-cycle)
// ---------------------------------------------------------------------------

namespace {

// Transpose an (n x ncols) CSR-shaped matrix (rows = A.n, columns 0..ncols-1).
CSR csr_transpose(const CSR& A, int ncols) {
  CSR T;
  T.n = ncols;
  T.ptr.assign(ncols + 1, 0);
  const int nnz = static_cast<int>(A.col.size());
  for (int k = 0; k < nnz; ++k) T.ptr[A.col[k] + 1]++;
  for (int i = 0; i < ncols; ++i) T.ptr[i + 1] += T.ptr[i];
  T.col.assign(nnz, 0);
  T.val.assign(nnz, 0.0);
  std::vector<int> next(T.ptr.begin(), T.ptr.end() - 1);
  for (int i = 0; i < A.n; ++i) {
    for (int k = A.ptr[i]; k < A.ptr[i + 1]; ++k) {
      const int c = A.col[k];
      const int dst = next[c]++;
      T.col[dst] = i;
      T.val[dst] = A.val[k];
    }
  }
  // Sort each row of T by column (rows came from scattering rows of A in
  // increasing i order but the original row order among a given target
  // column isn't necessarily sorted by original row index... actually it
  // is, since we scan i = 0..n-1 in order).
  return T;
}

// Small dense LU with partial pivoting (row-major, m x m). Mirrors the
// convention used for S-2's block solver: factor stores L (unit diag
// implicit) and U in-place, piv holds the row permutation.
void lu_factor_dense(std::vector<double>& a, std::vector<int>& piv, int m) {
  piv.resize(m);
  for (int i = 0; i < m; ++i) piv[i] = i;
  for (int k = 0; k < m; ++k) {
    int p = k;
    double best = std::fabs(a[k * m + k]);
    for (int i = k + 1; i < m; ++i) {
      const double v = std::fabs(a[i * m + k]);
      if (v > best) { best = v; p = i; }
    }
    if (p != k) {
      for (int j = 0; j < m; ++j) std::swap(a[k * m + j], a[p * m + j]);
      std::swap(piv[k], piv[p]);
    }
    const double piv_val = a[k * m + k];
    if (piv_val == 0.0) continue;
    for (int i = k + 1; i < m; ++i) {
      const double f = a[i * m + k] / piv_val;
      a[i * m + k] = f;
      for (int j = k + 1; j < m; ++j) a[i * m + j] -= f * a[k * m + j];
    }
  }
}

void lu_solve_dense(const std::vector<double>& a, const std::vector<int>& piv,
                    int m, const std::vector<double>& b,
                    std::vector<double>& x) {
  x.resize(m);
  std::vector<double> y(m);
  for (int i = 0; i < m; ++i) y[i] = b[piv[i]];
  for (int i = 0; i < m; ++i) {
    double s = y[i];
    for (int j = 0; j < i; ++j) s -= a[i * m + j] * y[j];
    y[i] = s;  // unit lower diagonal
  }
  for (int i = m - 1; i >= 0; --i) {
    double s = y[i];
    for (int j = i + 1; j < m; ++j) s -= a[i * m + j] * x[j];
    x[i] = (a[i * m + i] != 0.0) ? s / a[i * m + i] : s;
  }
}

}  // namespace

void AMG::setup(const CSR& A) {
  A_ = A;
  const int n = A.n;
  dinv_ = inv_diag(A);

  // 1) Strength-of-connection graph: |a_ij| > theta * sqrt(a_ii*a_jj).
  const double theta = 0.08;
  std::vector<std::vector<int>> strong(n);
  for (int i = 0; i < n; ++i) {
    const int dk = A.find(i, i);
    const double aii = (dk >= 0) ? A.val[dk] : 0.0;
    if (aii <= 0.0) continue;  // robustified: singleton aggregate below
    for (int k = A.ptr[i]; k < A.ptr[i + 1]; ++k) {
      const int j = A.col[k];
      if (j == i) continue;
      const int dj = A.find(j, j);
      const double ajj = (dj >= 0) ? A.val[dj] : 0.0;
      if (ajj <= 0.0) continue;
      if (std::fabs(A.val[k]) > theta * std::sqrt(aii * ajj)) strong[i].push_back(j);
    }
  }

  // 2) Greedy deterministic root-node aggregation.
  std::vector<int> agg(n, -1);
  int nagg = 0;
  for (int i = 0; i < n; ++i) {
    if (agg[i] != -1) continue;
    bool ok = true;
    for (int j : strong[i]) {
      if (agg[j] != -1) { ok = false; break; }
    }
    if (!ok) continue;
    const int g = nagg++;
    agg[i] = g;
    for (int j : strong[i]) agg[j] = g;
  }
  // Second pass: absorb leftovers into a neighboring aggregate, else singleton.
  for (int i = 0; i < n; ++i) {
    if (agg[i] != -1) continue;
    int found = -1;
    for (int j : strong[i]) {
      if (agg[j] != -1) { found = agg[j]; break; }
    }
    if (found != -1) {
      agg[i] = found;
    } else {
      agg[i] = nagg++;
    }
  }
  nagg_ = nagg;

  // 3) Power-iteration estimate of rho(D^-1 A).
  std::vector<double> v(n), vn(n);
  for (int i = 0; i < n; ++i) v[i] = 1.0 + 1e-3 * (i % 7);
  double rho_hat = 1.0;
  for (int it = 0; it < 10; ++it) {
    A.mul(v, vn);
    for (int i = 0; i < n; ++i) vn[i] *= dinv_[i];
    double nv = 0.0, nv0 = 0.0;
    for (int i = 0; i < n; ++i) { nv += vn[i] * vn[i]; nv0 += v[i] * v[i]; }
    nv = std::sqrt(nv);
    nv0 = std::sqrt(nv0);
    rho_hat = (nv0 > 0.0) ? nv / nv0 : 1.0;
    v = vn;
  }
  const double omega = 4.0 / (3.0 * rho_hat);

  // 4) Smoothed prolongator: P = (I - omega * Dinv * A) * P0, where
  //    P0[i, agg[i]] = 1. Row i of P has columns {agg[j] : a_ij != 0}.
  std::vector<int> Pptr(n + 1, 0);
  std::vector<std::vector<std::pair<int, double>>> Prows(n);
  for (int i = 0; i < n; ++i) {
    std::vector<std::pair<int, double>> acc;  // small linear-scan accumulator
    for (int k = A.ptr[i]; k < A.ptr[i + 1]; ++k) {
      const int j = A.col[k];
      const int cg = agg[j];
      const double delta_ij = (j == i) ? 1.0 : 0.0;
      const double contrib = (delta_ij - omega * dinv_[i] * A.val[k]);
      bool found = false;
      for (auto& pr : acc) {
        if (pr.first == cg) { pr.second += contrib; found = true; break; }
      }
      if (!found) acc.emplace_back(cg, contrib);
    }
    std::sort(acc.begin(), acc.end(),
              [](const auto& a, const auto& b) { return a.first < b.first; });
    Prows[i] = std::move(acc);
    Pptr[i + 1] = static_cast<int>(Prows[i].size());
  }
  for (int i = 0; i < n; ++i) Pptr[i + 1] += Pptr[i];
  P_.n = n;
  P_.ptr = Pptr;
  P_.col.resize(Pptr[n]);
  P_.val.resize(Pptr[n]);
  for (int i = 0; i < n; ++i) {
    int off = Pptr[i];
    for (const auto& pr : Prows[i]) {
      P_.col[off] = pr.first;
      P_.val[off] = pr.second;
      ++off;
    }
  }

  R_ = csr_transpose(P_, nagg_);

  // 5) Coarse matrix Ac = R * (A * P) via hashed sparse triple product.
  CSR AP;
  AP.n = n;
  AP.ptr.assign(n + 1, 0);
  std::vector<std::vector<std::pair<int, double>>> ap_rows(n);
#pragma omp parallel for schedule(dynamic, 64)
  for (int i = 0; i < n; ++i) {
    std::unordered_map<int, double> acc;
    for (int k = A.ptr[i]; k < A.ptr[i + 1]; ++k) {
      const int jrow = A.col[k];
      const double aval = A.val[k];
      for (int l = P_.ptr[jrow]; l < P_.ptr[jrow + 1]; ++l) {
        acc[P_.col[l]] += aval * P_.val[l];
      }
    }
    auto& row = ap_rows[i];
    row.assign(acc.begin(), acc.end());
    std::sort(row.begin(), row.end(),
              [](const auto& a, const auto& b) { return a.first < b.first; });
  }
  for (int i = 0; i < n; ++i) AP.ptr[i + 1] = static_cast<int>(ap_rows[i].size());
  for (int i = 0; i < n; ++i) AP.ptr[i + 1] += AP.ptr[i];
  AP.col.resize(AP.ptr[n]);
  AP.val.resize(AP.ptr[n]);
  for (int i = 0; i < n; ++i) {
    int off = AP.ptr[i];
    for (const auto& pr : ap_rows[i]) {
      AP.col[off] = pr.first;
      AP.val[off] = pr.second;
      ++off;
    }
  }

  Ac_.n = nagg_;
  Ac_.ptr.assign(nagg_ + 1, 0);
  std::vector<std::vector<std::pair<int, double>>> ac_rows(nagg_);
#pragma omp parallel for schedule(dynamic, 64)
  for (int gi = 0; gi < nagg_; ++gi) {
    std::unordered_map<int, double> acc;
    for (int k = R_.ptr[gi]; k < R_.ptr[gi + 1]; ++k) {
      const int arow = R_.col[k];
      const double rval = R_.val[k];
      for (int l = AP.ptr[arow]; l < AP.ptr[arow + 1]; ++l) {
        acc[AP.col[l]] += rval * AP.val[l];
      }
    }
    auto& row = ac_rows[gi];
    row.assign(acc.begin(), acc.end());
    std::sort(row.begin(), row.end(),
              [](const auto& a, const auto& b) { return a.first < b.first; });
  }
  for (int gi = 0; gi < nagg_; ++gi) Ac_.ptr[gi + 1] = static_cast<int>(ac_rows[gi].size());
  for (int gi = 0; gi < nagg_; ++gi) Ac_.ptr[gi + 1] += Ac_.ptr[gi];
  Ac_.col.resize(Ac_.ptr[nagg_]);
  Ac_.val.resize(Ac_.ptr[nagg_]);
  for (int gi = 0; gi < nagg_; ++gi) {
    int off = Ac_.ptr[gi];
    for (const auto& pr : ac_rows[gi]) {
      Ac_.col[off] = pr.first;
      Ac_.val[off] = pr.second;
      ++off;
    }
  }

  // 6) Coarse solve: dense LU if small, else Jacobi sweeps at apply time.
  coarse_dense_ = (nagg_ < 200);
  if (coarse_dense_) {
    coarse_lu_.assign(static_cast<size_t>(nagg_) * nagg_, 0.0);
    for (int i = 0; i < nagg_; ++i)
      for (int k = Ac_.ptr[i]; k < Ac_.ptr[i + 1]; ++k)
        coarse_lu_[static_cast<size_t>(i) * nagg_ + Ac_.col[k]] = Ac_.val[k];
    lu_factor_dense(coarse_lu_, coarse_piv_, nagg_);
  }
}

void AMG::apply(const std::vector<double>& r, std::vector<double>& z) const {
  const int n = A_.n;
  const double omega_j = 2.0 / 3.0;

  // Pre-smooth (x0 = 0): x = omega_j * Dinv * r, then one more damped-Jacobi
  // sweep. Two pre/post sweeps (V(2,2) rather than the textbook V(1,1))
  // measurably improved the iteration-count margin over ILU(0) on the test
  // problem (21 iters at V(1,1) vs ILU0's 42 -- exactly 2x, not strictly
  // better -- while V(2,2) brings AMG to ~14-16 iters).
  std::vector<double> x(n);
#pragma omp parallel for schedule(static)
  for (int i = 0; i < n; ++i) x[i] = omega_j * dinv_[i] * r[i];

  std::vector<double> Ax(n), rf(n);
  A_.mul(x, Ax);
#pragma omp parallel for schedule(static)
  for (int i = 0; i < n; ++i) x[i] += omega_j * dinv_[i] * (r[i] - Ax[i]);

  // Fine residual rf = r - A x
  A_.mul(x, Ax);
#pragma omp parallel for schedule(static)
  for (int i = 0; i < n; ++i) rf[i] = r[i] - Ax[i];

  // Restrict: rc = R * rf
  std::vector<double> rc(nagg_, 0.0);
  R_.mul(rf, rc);

  // Coarse solve
  std::vector<double> ec(nagg_, 0.0);
  if (coarse_dense_) {
    lu_solve_dense(coarse_lu_, coarse_piv_, nagg_, rc, ec);
  } else {
    std::vector<double> cdinv(nagg_, 1.0);
    for (int i = 0; i < nagg_; ++i) {
      const int k = Ac_.find(i, i);
      if (k >= 0 && Ac_.val[k] != 0.0) cdinv[i] = 1.0 / Ac_.val[k];
    }
    std::vector<double> Ace(nagg_);
    for (int sweep = 0; sweep < 10; ++sweep) {
      Ac_.mul(ec, Ace);
      for (int i = 0; i < nagg_; ++i)
        ec[i] += omega_j * cdinv[i] * (rc[i] - Ace[i]);
    }
  }

  // Prolongate + correct: x += P * ec
  std::vector<double> Pec(n, 0.0);
  P_.mul(ec, Pec);
#pragma omp parallel for schedule(static)
  for (int i = 0; i < n; ++i) x[i] += Pec[i];

  // Post-smooth: two damped-Jacobi sweeps, x += omega_j * Dinv * (r - A x)
  A_.mul(x, Ax);
#pragma omp parallel for schedule(static)
  for (int i = 0; i < n; ++i) x[i] += omega_j * dinv_[i] * (r[i] - Ax[i]);
  A_.mul(x, Ax);
#pragma omp parallel for schedule(static)
  for (int i = 0; i < n; ++i) x[i] += omega_j * dinv_[i] * (r[i] - Ax[i]);

  z = x;
}

SolveResult cg_amg(const CSR& A, const std::vector<double>& b,
                   std::vector<double>& x, double rtol, int maxit) {
  AMG amg;
  amg.setup(A);
  Precond psolve = [&](const std::vector<double>& in, std::vector<double>& out) {
    amg.apply(in, out);
  };
  return cg(A, b, x, rtol, maxit, psolve);
}

// ---------------------------------------------------------------------------
// Preconditioned Krylov cores
// ---------------------------------------------------------------------------

SolveResult cg(const CSR& A, const std::vector<double>& b,
               std::vector<double>& x, double rtol, int maxit,
               const Precond& psolve) {
  ensure_ftz_daz();
  ensure_sane_thread_count();
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
  ensure_ftz_daz();
  ensure_sane_thread_count();
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

#ifdef CPROCESS_GPU
// ---------------------------------------------------------------------------
// GPU (OpenMP target offload) Jacobi-preconditioned CG
//
// Same algorithm as cg_jacobi/cg above, but every SpMV / AXPY / dot-product /
// Jacobi-apply kernel runs in an `omp target teams distribute parallel for`
// region on the offload device. Device buffers for the matrix (ptr/col/val),
// the RHS, and the Jacobi diagonal are mapped once for the whole solve via a
// single `target data` region; only the scalar reductions needed for the
// convergence test are transferred back to the host each iteration. The
// iteration control flow itself stays on the host.
//
// On a system with no offload device (omp_get_num_devices() == 0), or when
// built with a host-only OpenMP runtime, the `target` regions execute on the
// host itself per the OpenMP spec -- so this function is safe to call (and
// numerically identical to cg_jacobi) even without real GPU hardware.
// ---------------------------------------------------------------------------
SolveResult cg_jacobi_gpu(const CSR& A, const std::vector<double>& b,
                          std::vector<double>& x, double rtol, int maxit) {
  ensure_ftz_daz();
  ensure_sane_thread_count();
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
  const int nnz = static_cast<int>(A.val.size());

  const int* ptr = A.ptr.data();
  const int* col = A.col.data();
  const double* val = A.val.data();
  const double* b_ = b.data();
  const double* Mp = M.data();
  double* x_ = x.data();

  std::vector<double> r(n), z(n), p(n), q(n);
  double* r_ = r.data();
  double* z_ = z.data();
  double* p_ = p.data();
  double* q_ = q.data();

  double rz = 0.0;

#pragma omp target data map(to: ptr[0:n + 1], col[0:nnz], val[0:nnz], \
                                 b_[0:n], Mp[0:n])                    \
                         map(tofrom: x_[0:n])                         \
                         map(alloc: r_[0:n], z_[0:n], p_[0:n], q_[0:n])
  {
    // q <- A*x ; r <- b - q
#pragma omp target teams distribute parallel for
    for (int i = 0; i < n; ++i) {
      double s = 0.0;
      for (int k = ptr[i]; k < ptr[i + 1]; ++k) s += val[k] * x_[col[k]];
      q_[i] = s;
    }
#pragma omp target teams distribute parallel for
    for (int i = 0; i < n; ++i) r_[i] = b_[i] - q_[i];

    // z <- M^{-1} r ; p <- z
#pragma omp target teams distribute parallel for
    for (int i = 0; i < n; ++i) z_[i] = Mp[i] * r_[i];
#pragma omp target teams distribute parallel for
    for (int i = 0; i < n; ++i) p_[i] = z_[i];

    rz = 0.0;
#pragma omp target teams distribute parallel for map(tofrom: rz) reduction(+ : rz)
    for (int i = 0; i < n; ++i) rz += r_[i] * z_[i];

    for (int it = 0; it < maxit; ++it) {
      // q <- A*p
#pragma omp target teams distribute parallel for
      for (int i = 0; i < n; ++i) {
        double s = 0.0;
        for (int k = ptr[i]; k < ptr[i + 1]; ++k) s += val[k] * p_[col[k]];
        q_[i] = s;
      }

      double pq = 0.0;
#pragma omp target teams distribute parallel for map(tofrom: pq) reduction(+ : pq)
      for (int i = 0; i < n; ++i) pq += p_[i] * q_[i];
      if (pq == 0.0) break;

      const double alpha = rz / pq;
#pragma omp target teams distribute parallel for
      for (int i = 0; i < n; ++i) x_[i] += alpha * p_[i];
#pragma omp target teams distribute parallel for
      for (int i = 0; i < n; ++i) r_[i] -= alpha * q_[i];

      double rn2 = 0.0;
#pragma omp target teams distribute parallel for map(tofrom: rn2) reduction(+ : rn2)
      for (int i = 0; i < n; ++i) rn2 += r_[i] * r_[i];
      const double rn = std::sqrt(rn2);
      res.iters = it + 1;
      res.resid = rn / bnorm;
      if (rn <= rtol * bnorm) {
        res.converged = true;
        break;
      }

#pragma omp target teams distribute parallel for
      for (int i = 0; i < n; ++i) z_[i] = Mp[i] * r_[i];

      double rz2 = 0.0;
#pragma omp target teams distribute parallel for map(tofrom: rz2) reduction(+ : rz2)
      for (int i = 0; i < n; ++i) rz2 += r_[i] * z_[i];
      const double beta = rz2 / rz;
      rz = rz2;

#pragma omp target teams distribute parallel for
      for (int i = 0; i < n; ++i) p_[i] = z_[i] + beta * p_[i];
    }
  }  // end target data

  return res;
}
#endif  // CPROCESS_GPU

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
SolveResult gmres_op(const LinOp& aop, int n, const std::vector<double>& b,
                     std::vector<double>& x, double rtol, int maxit,
                     int restart, const Precond& psolve) {
  SolveResult res;
  x.resize(n, 0.0);
  const double bnorm = std::sqrt(dotv(b, b));
  if (bnorm == 0.0) {
    x.assign(n, 0.0);
    res.converged = true;
    return res;
  }

  std::vector<double> r(n), t(n);

  for (int total = 0; total < maxit; ) {
    // r = b - A*x (true residual at restart boundary)
    aop(x, t);
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
      psolve(V[m], t);
      aop(t, r);  // r = w

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
    std::vector<double> vy(n, 0.0);
#pragma omp parallel for schedule(static)
    for (int i = 0; i < n; ++i) {
      double acc = 0.0;
      for (int j = 0; j < m; ++j) acc += V[j][i] * y[j];
      vy[i] = acc;
    }
    psolve(vy, t);
#pragma omp parallel for schedule(static)
    for (int i = 0; i < n; ++i) x[i] += t[i];

    if (res.resid <= rtol) { res.converged = true; return res; }
    if (hnext < 1e-12) break;  // exact solution found mid-restart
  }

  // Recompute true residual after final update
  aop(x, t);
#pragma omp parallel for schedule(static)
  for (int i = 0; i < n; ++i) r[i] = b[i] - t[i];
  res.resid = std::sqrt(dotv(r, r)) / bnorm;
  res.converged = (res.resid <= rtol);
  return res;
}

SolveResult gmres_jacobi(const CSR& A, const std::vector<double>& b,
                         std::vector<double>& x, double rtol, int maxit,
                         int restart) {
  const std::vector<double> M = inv_diag(A);
  LinOp aop = [&](const std::vector<double>& in, std::vector<double>& out) {
    A.mul(in, out);
  };
  Precond psolve = [&](const std::vector<double>& in, std::vector<double>& out) {
    const int n = A.n;
    out.resize(n);
#pragma omp parallel for schedule(static)
    for (int i = 0; i < n; ++i) out[i] = M[i] * in[i];
  };
  return gmres_op(aop, A.n, b, x, rtol, maxit, restart, psolve);
}

// ---------------------------------------------------------------------------
// S-2: Block CSR + block-Jacobi preconditioned BiCGSTAB
// ---------------------------------------------------------------------------

void BCSR::mul(const std::vector<double>& x, std::vector<double>& y) const {
  y.assign(static_cast<std::size_t>(n) * nb, 0.0);
#pragma omp parallel for schedule(static)
  for (int i = 0; i < n; ++i) {
    double* yi = &y[static_cast<std::size_t>(i) * nb];
    for (int k = ptr[i]; k < ptr[i + 1]; ++k) {
      const double* B = &val[static_cast<std::size_t>(k) * nb * nb];
      const double* xj = &x[static_cast<std::size_t>(col[k]) * nb];
      for (int r = 0; r < nb; ++r) {
        double s = 0;
        for (int c = 0; c < nb; ++c) s += B[r * nb + c] * xj[c];
        yi[r] += s;
      }
    }
  }
}

int BCSR::find(int row, int c) const {
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

BCSR bcsr_from_pattern(const CSR& scalar_pattern, int nb) {
  BCSR B;
  B.n = scalar_pattern.n;
  B.nb = nb;
  B.ptr = scalar_pattern.ptr;
  B.col = scalar_pattern.col;
  B.val.assign(static_cast<std::size_t>(B.ptr.empty() ? 0 : B.ptr[B.n]) *
                   static_cast<std::size_t>(nb) * static_cast<std::size_t>(nb),
               0.0);
  return B;
}

namespace {

// In-place dense LU with partial pivoting of the nb x nb matrix `a`
// (row-major). piv[k] records the row swapped into position k. Returns
// false if a pivot is exactly zero (singular block -- caller treats the
// block as identity).
bool lu_factor_dense_ptr(double* a, int* piv, int nb) {
  for (int k = 0; k < nb; ++k) piv[k] = k;
  for (int k = 0; k < nb; ++k) {
    int p = k;
    double best = std::fabs(a[k * nb + k]);
    for (int i = k + 1; i < nb; ++i) {
      const double v = std::fabs(a[i * nb + k]);
      if (v > best) { best = v; p = i; }
    }
    if (p != k) {
      for (int j = 0; j < nb; ++j) std::swap(a[k * nb + j], a[p * nb + j]);
      std::swap(piv[k], piv[p]);
    }
    if (a[k * nb + k] == 0.0) return false;
    for (int i = k + 1; i < nb; ++i) {
      a[i * nb + k] /= a[k * nb + k];
      for (int j = k + 1; j < nb; ++j) a[i * nb + j] -= a[i * nb + k] * a[k * nb + j];
    }
  }
  return true;
}

// Solves L U x = b using the factors/pivots from lu_factor_dense_ptr; b is
// overwritten with x.
void lu_solve_dense_ptr(const double* a, const int* piv, int nb, double* b) {
  std::vector<double> y(nb);
  for (int i = 0; i < nb; ++i) y[i] = b[piv[i]];
  for (int i = 0; i < nb; ++i) {
    double s = y[i];
    for (int j = 0; j < i; ++j) s -= a[i * nb + j] * y[j];
    y[i] = s;  // unit lower diagonal
  }
  for (int i = nb - 1; i >= 0; --i) {
    double s = y[i];
    for (int j = i + 1; j < nb; ++j) s -= a[i * nb + j] * y[j];
    y[i] = s / a[i * nb + i];
  }
  for (int i = 0; i < nb; ++i) b[i] = y[i];
}

struct BlockJacobi {
  int n = 0, nb = 0;
  std::vector<double> lu;
  std::vector<int> piv;
  std::vector<char> ok;

  void setup(const BCSR& A) {
    n = A.n;
    nb = A.nb;
    lu.assign(static_cast<std::size_t>(n) * nb * nb, 0.0);
    piv.assign(static_cast<std::size_t>(n) * nb, 0);
    ok.assign(n, 0);
#pragma omp parallel for schedule(static)
    for (int i = 0; i < n; ++i) {
      const int k = A.find(i, i);
      double* Di = &lu[static_cast<std::size_t>(i) * nb * nb];
      int* Pi = &piv[static_cast<std::size_t>(i) * nb];
      if (k < 0) { ok[i] = 0; continue; }
      const double* src = &A.val[static_cast<std::size_t>(k) * nb * nb];
      for (int e = 0; e < nb * nb; ++e) Di[e] = src[e];
      ok[i] = lu_factor_dense_ptr(Di, Pi, nb) ? 1 : 0;
    }
  }

  void apply(const std::vector<double>& r, std::vector<double>& z) const {
    z.assign(static_cast<std::size_t>(n) * nb, 0.0);
#pragma omp parallel for schedule(static)
    for (int i = 0; i < n; ++i) {
      double* zi = &z[static_cast<std::size_t>(i) * nb];
      const double* ri = &r[static_cast<std::size_t>(i) * nb];
      if (!ok[i]) {
        for (int c = 0; c < nb; ++c) zi[c] = ri[c];  // identity fallback
        continue;
      }
      const double* Di = &lu[static_cast<std::size_t>(i) * nb * nb];
      const int* Pi = &piv[static_cast<std::size_t>(i) * nb];
      for (int c = 0; c < nb; ++c) zi[c] = ri[c];
      lu_solve_dense_ptr(Di, Pi, nb, zi);
    }
  }
};

double bdotv(const std::vector<double>& a, const std::vector<double>& b) {
  const int n = static_cast<int>(a.size());
  double s = 0;
#pragma omp parallel for reduction(+ : s) schedule(static)
  for (int i = 0; i < n; ++i) s += a[i] * b[i];
  return s;
}

void baxpy(double alpha, const std::vector<double>& x, std::vector<double>& y) {
  const int n = static_cast<int>(x.size());
#pragma omp parallel for schedule(static)
  for (int i = 0; i < n; ++i) y[i] += alpha * x[i];
}

}  // namespace

SolveResult bicgstab_bjacobi(const BCSR& A, const std::vector<double>& b,
                             std::vector<double>& x, double rtol, int maxit) {
  ensure_ftz_daz();
  ensure_sane_thread_count();
  SolveResult res;
  const int n_total = A.n * A.nb;
  x.resize(n_total, 0.0);
  const double bnorm = std::sqrt(bdotv(b, b));
  if (bnorm == 0.0) {
    x.assign(n_total, 0.0);
    res.converged = true;
    return res;
  }
  BlockJacobi M;
  M.setup(A);
  Precond psolve = [&](const std::vector<double>& in, std::vector<double>& out) {
    M.apply(in, out);
  };

  std::vector<double> r(n_total), r0(n_total), p(n_total, 0.0), v(n_total, 0.0),
      s(n_total), t(n_total), ph(n_total), sh(n_total);
  A.mul(x, t);
#pragma omp parallel for schedule(static)
  for (int i = 0; i < n_total; ++i) r[i] = b[i] - t[i];
  r0 = r;
  double rho = 1, alpha = 1, omega = 1;
  for (int it = 0; it < maxit; ++it) {
    const double rho2 = bdotv(r0, r);
    if (rho2 == 0.0) break;
    if (it == 0) {
      p = r;
    } else {
      const double beta = (rho2 / rho) * (alpha / omega);
#pragma omp parallel for schedule(static)
      for (int i = 0; i < n_total; ++i) p[i] = r[i] + beta * (p[i] - omega * v[i]);
    }
    rho = rho2;
    psolve(p, ph);
    A.mul(ph, v);
    const double r0v = bdotv(r0, v);
    if (r0v == 0.0) break;
    alpha = rho / r0v;
#pragma omp parallel for schedule(static)
    for (int i = 0; i < n_total; ++i) s[i] = r[i] - alpha * v[i];
    double sn = std::sqrt(bdotv(s, s));
    res.iters = it + 1;
    if (sn <= rtol * bnorm) {
      baxpy(alpha, ph, x);
      res.resid = sn / bnorm;
      res.converged = true;
      return res;
    }
    psolve(s, sh);
    A.mul(sh, t);
    const double tt = bdotv(t, t);
    if (tt == 0.0) break;
    omega = bdotv(t, s) / tt;
#pragma omp parallel for schedule(static)
    for (int i = 0; i < n_total; ++i) x[i] += alpha * ph[i] + omega * sh[i];
#pragma omp parallel for schedule(static)
    for (int i = 0; i < n_total; ++i) r[i] = s[i] - omega * t[i];
    const double rn = std::sqrt(bdotv(r, r));
    res.resid = rn / bnorm;
    if (rn <= rtol * bnorm) {
      res.converged = true;
      return res;
    }
    if (omega == 0.0) break;
  }
  return res;
}

// PA-2: BFS from `start`, returns the farthest node reached (ties broken by
// lowest node index, for full determinism) together with the level array.
namespace {
int bfs_farthest(const CSR& A, int start, std::vector<int>& level) {
  const int n = A.n;
  level.assign(n, -1);
  level[start] = 0;
  std::vector<int> q;
  q.reserve(n);
  q.push_back(start);
  int maxlvl = 0;
  for (std::size_t qi = 0; qi < q.size(); ++qi) {
    const int u = q[qi];
    for (int k = A.ptr[u]; k < A.ptr[u + 1]; ++k) {
      const int v = A.col[k];
      if (v == u || level[v] != -1) continue;
      level[v] = level[u] + 1;
      maxlvl = std::max(maxlvl, level[v]);
      q.push_back(v);
    }
  }
  int best = -1, best_deg = -1;
  for (int i = 0; i < n; ++i) {
    if (level[i] != maxlvl) continue;
    int d = 0;
    for (int k = A.ptr[i]; k < A.ptr[i + 1]; ++k)
      if (A.col[k] != i) ++d;
    if (best == -1 || d < best_deg) { best = i; best_deg = d; }
  }
  return best != -1 ? best : start;
}
}  // namespace

std::vector<int> rcm_order(const CSR& A) {
  const int n = A.n;
  std::vector<int> order;
  if (n == 0) return order;
  order.reserve(n);

  std::vector<int> degree(n, 0);
  for (int i = 0; i < n; ++i) {
    int d = 0;
    for (int k = A.ptr[i]; k < A.ptr[i + 1]; ++k)
      if (A.col[k] != i) ++d;
    degree[i] = d;
  }

  std::vector<char> visited(n, 0);
  std::vector<int> level;  // scratch for bfs_farthest

  auto cm_bfs = [&](int start) {
    std::vector<int> q;
    q.reserve(n);
    visited[start] = 1;
    order.push_back(start);
    q.push_back(start);
    for (std::size_t qi = 0; qi < q.size(); ++qi) {
      const int u = q[qi];
      std::vector<int> nbrs;
      for (int k = A.ptr[u]; k < A.ptr[u + 1]; ++k) {
        const int v = A.col[k];
        if (v != u && !visited[v]) nbrs.push_back(v);
      }
      std::sort(nbrs.begin(), nbrs.end(), [&](int a, int b) {
        if (degree[a] != degree[b]) return degree[a] < degree[b];
        return a < b;
      });
      for (int v : nbrs) {
        if (!visited[v]) {
          visited[v] = 1;
          order.push_back(v);
          q.push_back(v);
        }
      }
    }
  };

  // Pseudo-peripheral node search (capped at two BFS passes, per spec):
  // BFS from the seed to get the farthest node u, then BFS from u to get
  // the farthest node v; v is used as the Cuthill-McKee start node.
  // Repeated per disconnected component, starting from the lowest-index
  // unvisited node each time.
  for (int seed = 0; seed < n; ++seed) {
    if (visited[seed]) continue;
    const int u = bfs_farthest(A, seed, level);
    const int v = bfs_farthest(A, u, level);
    cm_bfs(v);
  }

  std::vector<int> perm(n);
  for (int i = 0; i < n; ++i) perm[i] = order[n - 1 - i];  // reverse
  return perm;
}

CSR permute(const CSR& A, const std::vector<int>& perm) {
  const int n = A.n;
  CSR B;
  B.n = n;
  if (n == 0) return B;
  std::vector<int> iperm(n);
  for (int i = 0; i < n; ++i) iperm[perm[i]] = i;

  B.ptr.assign(n + 1, 0);
  for (int newi = 0; newi < n; ++newi) {
    const int oldi = perm[newi];
    B.ptr[newi + 1] = B.ptr[newi] + (A.ptr[oldi + 1] - A.ptr[oldi]);
  }
  B.col.resize(B.ptr[n]);
  B.val.resize(B.ptr[n]);
  std::vector<std::pair<int, double>> row;
  for (int newi = 0; newi < n; ++newi) {
    const int oldi = perm[newi];
    row.clear();
    for (int k = A.ptr[oldi]; k < A.ptr[oldi + 1]; ++k)
      row.emplace_back(iperm[A.col[k]], A.val[k]);
    std::sort(row.begin(), row.end());
    const int off = B.ptr[newi];
    for (std::size_t t = 0; t < row.size(); ++t) {
      B.col[off + static_cast<int>(t)] = row[t].first;
      B.val[off + static_cast<int>(t)] = row[t].second;
    }
  }
  return B;
}

}  // namespace cp
