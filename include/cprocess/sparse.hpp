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

// Jacobi-preconditioned BiCGSTAB (fallback for non-symmetric systems).
SolveResult bicgstab_jacobi(const CSR& A, const std::vector<double>& b,
                            std::vector<double>& x, double rtol, int maxit);

// Jacobi-preconditioned restarted GMRES(restart).
// Arnoldi + Givens rotations, right-preconditioned (M^{-1} = Jacobi).
SolveResult gmres_jacobi(const CSR& A, const std::vector<double>& b,
                         std::vector<double>& x, double rtol, int maxit,
                         int restart = 30);

// ILU(0)-preconditioned variants (parallel triangular solves). The factor is
// built internally from A on every call; pass a prebuilt one via the *_pc forms
// to reuse the factorization across right-hand sides.
SolveResult cg_ilu0(const CSR& A, const std::vector<double>& b,
                    std::vector<double>& x, double rtol, int maxit);
SolveResult bicgstab_ilu0(const CSR& A, const std::vector<double>& b,
                          std::vector<double>& x, double rtol, int maxit);

// Generic cores taking an explicit preconditioner (M^{-1} apply).
SolveResult cg(const CSR& A, const std::vector<double>& b,
               std::vector<double>& x, double rtol, int maxit,
               const Precond& psolve);
SolveResult bicgstab(const CSR& A, const std::vector<double>& b,
                     std::vector<double>& x, double rtol, int maxit,
                     const Precond& psolve);

}  // namespace cp
