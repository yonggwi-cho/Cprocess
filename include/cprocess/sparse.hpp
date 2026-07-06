#pragma once
#include <functional>
#include <vector>

namespace cp {

// Compressed sparse row matrix (square).
struct CSR {
  int n = 0;
  std::vector<int> ptr;     // size n+1
  std::vector<int> col;     // size nnz, column indices sorted within a row
  std::vector<double> val;  // size nnz

  void mul(const std::vector<double>& x, std::vector<double>& y) const;
  int find(int row, int c) const;  // value slot of (row, c), -1 if absent
};

// PA-2: Reverse Cuthill-McKee ordering of the (symmetric-pattern) matrix A,
// to shrink matrix bandwidth for better SpMV cache locality. Returns
// perm with perm[new] = old. Deterministic: pseudo-peripheral node search
// (BFS from node 0, then BFS from the farthest node found) picks the CM
// start node; neighbors are queued in (degree asc, index asc) order;
// disconnected components are each processed in turn, starting from the
// lowest-index unvisited node.
std::vector<int> rcm_order(const CSR& A);

// Symmetric permutation: B(new_i, new_j) = A(perm[new_i], perm[new_j]).
// Row columns of the result are sorted by column index.
CSR permute(const CSR& A, const std::vector<int>& perm);

struct SolveResult {
  bool converged = false;
  int iters = 0;
  double resid = 0;  // final ||b - A x|| / ||b||
};

// Preconditioner application: y <- M^{-1} x. Must not alias x and y.
using Precond = std::function<void(const std::vector<double>& x,
                                   std::vector<double>& y)>;

// ILU(0) factorization with level-scheduled parallel triangular solves.
//
// The incomplete factor L\U is stored in the sparsity pattern of A (no
// fill-in). Forward/backward substitution is parallelized across "levels":
// rows whose off-diagonal dependencies are all in earlier levels can be
// solved simultaneously. Level construction is done once at factor() time.
struct ILU0 {
  CSR lu;                                 // combined L (unit lower) and U
  std::vector<int> diag;                  // diagonal slot per row in `lu`
  std::vector<std::vector<int>> lvl_lo;   // level sets for the lower solve
  std::vector<std::vector<int>> lvl_up;   // level sets for the upper solve

  void factor(const CSR& A);              // compute the incomplete factors
  void apply(const std::vector<double>& x, std::vector<double>& y) const;
};

// Jacobi-preconditioned conjugate gradient (A must be SPD).
SolveResult cg_jacobi(const CSR& A, const std::vector<double>& b,
                      std::vector<double>& x, double rtol, int maxit);

#ifdef CPROCESS_GPU
// GPU (OpenMP target offload) Jacobi-preconditioned CG. Numerically
// equivalent to cg_jacobi (same algorithm), but SpMV/AXPY/dot/Jacobi
// kernels execute on the offload device; the iteration itself stays on
// the host (only reduction scalars are transferred back per iteration).
// Falls back to host execution when no offload device is present
// (omp_get_num_devices() == 0), so it is safe to call in CPU-only
// development environments as long as CPROCESS_GPU was enabled at
// configure time. Not declared at all when CPROCESS_GPU is undefined.
SolveResult cg_jacobi_gpu(const CSR& A, const std::vector<double>& b,
                          std::vector<double>& x, double rtol, int maxit);
#endif

// Jacobi-preconditioned BiCGSTAB (fallback for non-symmetric systems).
SolveResult bicgstab_jacobi(const CSR& A, const std::vector<double>& b,
                            std::vector<double>& x, double rtol, int maxit);

// Jacobi-preconditioned restarted GMRES(restart).
// Arnoldi + Givens rotations, right-preconditioned (M^{-1} = Jacobi).
SolveResult gmres_jacobi(const CSR& A, const std::vector<double>& b,
                         std::vector<double>& x, double rtol, int maxit,
                         int restart = 30);

// Generic operator form of restarted GMRES: right-preconditioned, works with
// any linear operator `aop` (y <- A x) and preconditioner `psolve`
// (y <- M^{-1} x) of dimension `n` -- no CSR required (used by S-3's
// Jacobian-free Newton-Krylov, where the "matrix" is a finite-difference
// directional derivative). `gmres_jacobi` above is a thin wrapper around
// this with aop = A.mul and psolve = Jacobi(inv_diag(A)).
using LinOp = std::function<void(const std::vector<double>&,
                                 std::vector<double>&)>;
SolveResult gmres_op(const LinOp& aop, int n, const std::vector<double>& b,
                     std::vector<double>& x, double rtol, int maxit,
                     int restart, const Precond& psolve);

// ILU(0)-preconditioned variants (parallel triangular solves). The factor is
// built internally from A on every call; pass a prebuilt one via the *_pc forms
// to reuse the factorization across right-hand sides.
SolveResult cg_ilu0(const CSR& A, const std::vector<double>& b,
                    std::vector<double>& x, double rtol, int maxit);
SolveResult bicgstab_ilu0(const CSR& A, const std::vector<double>& b,
                          std::vector<double>& x, double rtol, int maxit);

// Smoothed-aggregation AMG, two-level V-cycle. Usable as a Precond for
// the existing cg/bicgstab cores.
struct AMG {
  void setup(const CSR& A);   // build aggregates, P, R=P^T, Ac
  void apply(const std::vector<double>& r, std::vector<double>& z) const;
  int n_aggregates() const { return nagg_; }

 private:
  CSR A_;          // fine matrix (copy)
  CSR P_, R_;      // prolongation (n x nagg) and restriction (= P^T)
  CSR Ac_;         // coarse matrix R A P (nagg x nagg)
  std::vector<double> dinv_;       // fine 1/a_ii for the Jacobi smoother
  std::vector<double> coarse_lu_;  // dense LU of Ac when nagg < 200
  std::vector<int> coarse_piv_;
  bool coarse_dense_ = false;
  int nagg_ = 0;
};

SolveResult cg_amg(const CSR& A, const std::vector<double>& b,
                   std::vector<double>& x, double rtol, int maxit);

// Generic cores taking an explicit preconditioner (M^{-1} apply).
SolveResult cg(const CSR& A, const std::vector<double>& b,
               std::vector<double>& x, double rtol, int maxit,
               const Precond& psolve);
SolveResult bicgstab(const CSR& A, const std::vector<double>& b,
                     std::vector<double>& x, double rtol, int maxit,
                     const Precond& psolve);

// Block CSR: square block matrix of n block-rows; each stored entry is a
// dense nb x nb block (row-major). val.size() == ptr[n] * nb * nb.
struct BCSR {
  int n = 0;    // block rows
  int nb = 0;   // block size (number of coupled species)
  std::vector<int> ptr;      // size n+1 (block entries)
  std::vector<int> col;      // size nnz_blocks, sorted within a row
  std::vector<double> val;   // nnz_blocks * nb*nb, block k at val[k*nb*nb]

  // y <- A x, with x/y of length n*nb (cell-major: x[i*nb + s]).
  void mul(const std::vector<double>& x, std::vector<double>& y) const;
  int find(int row, int c) const;  // block slot index, -1 if absent
};

// Builds a BCSR sharing the sparsity pattern of a scalar CSR (val zeroed).
BCSR bcsr_from_pattern(const CSR& scalar_pattern, int nb);

// Block-Jacobi preconditioned BiCGSTAB on a BCSR system.
SolveResult bicgstab_bjacobi(const BCSR& A, const std::vector<double>& b,
                             std::vector<double>& x, double rtol, int maxit);

}  // namespace cp
