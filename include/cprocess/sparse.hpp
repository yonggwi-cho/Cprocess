#pragma once
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

}  // namespace cp
