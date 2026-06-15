#pragma once
#include <array>
#include <iosfwd>
#include <string>
#include <vector>

#include "materials.hpp"
#include "mesh.hpp"
#include "sparse.hpp"

namespace cp {

// Dirichlet (fixed surface concentration) condition on a boundary patch,
// e.g. a pre-deposition source. Default everywhere else: zero flux.
struct DirichletBC {
  std::string species;  // dopant symbol, e.g. "P"
  int patch = -1;
  double conc = 0;      // cm^-3
};

struct DiffuseOpts {
  double temp = 1273.15;   // K
  double time = 60;        // s
  double dt = 0;           // s; <=0 selects time/50
  bool field_enh = true;   // electric-field drift enhancement (factor <= 2)
  bool nonortho = true;    // deferred non-orthogonal flux correction
  int max_picard = 8;
  double picard_tol = 1e-4;
  double lin_rtol = 1e-10;
  int lin_maxit = 2000;
  int verbosity = 1;       // 0 silent, 1 per-step lines
};

struct SpeciesField {
  const Dopant* dopant = nullptr;
  std::vector<double>* conc = nullptr;  // per-cell concentration, cm^-3
};

// Cell-centered finite-volume dopant diffusion on an unstructured tet mesh.
//
//   dC/dt = div( D(T, n/ni) grad C )
//
// Discretization: two-point flux with harmonic face diffusivity plus a
// deferred (explicit) non-orthogonal correction from weighted least-squares
// cell gradients (exact for linear fields on arbitrary tet meshes);
// backward-Euler time stepping; Picard iteration over the concentration-
// dependent diffusivity which couples all species through n/ni;
// Jacobi-preconditioned CG for the linear systems.
//
// Cells with solve_mask == 0 (non-silicon regions) are frozen and their
// interfaces are treated as zero-flux walls.
class DiffusionSolver {
 public:
  DiffusionSolver(const Mesh& mesh, std::vector<char> solve_mask,
                  std::ostream* log = nullptr);

  void run(std::vector<SpeciesField>& fields,
           const std::vector<DirichletBC>& bcs, const DiffuseOpts& opts);

 private:
  enum FaceKind : char {
    kInactive = 0,   // both sides masked
    kInternal = 1,   // flux between two active cells
    kBoundOwner = 2, // only owner active (boundary or material interface)
    kBoundNeigh = 3, // only neigh active (material interface)
  };
  struct FGeom {
    char kind = kInactive;
    double g = 0;             // two-point geometric transmissibility |S|^2/(S.d)
    double wP = 0;            // face interpolation weight of the owner cell
    double delP = 0, delN = 0;  // centroid-to-face distances along the normal
    double gb = 0;            // boundary transmissibility |S|^2/(S.(cf-cP))
    Vec3 k;                   // non-orthogonal residual vector S - g*d
  };

  void build();
  void gradients(const std::vector<double>& c, const std::vector<double>& bcface,
                 std::vector<Vec3>& grad) const;

  const Mesh& mesh_;
  std::vector<char> mask_;
  std::ostream* log_;
  std::vector<FGeom> fg_;
  CSR A_;
  std::vector<int> diag_;                  // diagonal slot per cell
  std::vector<std::array<int, 2>> fslot_;  // (P,N) and (N,P) slots per face
  int clamped_faces_ = 0;                  // faces with poor orthogonality
};

}  // namespace cp
