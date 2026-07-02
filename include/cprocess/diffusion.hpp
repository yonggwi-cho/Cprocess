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
  bool activation = true;  // clamp charge neutrality at solid solubility (P1-3)
};

struct SpeciesField {
  const Dopant* dopant = nullptr;
  std::vector<double>* conc = nullptr;  // per-cell concentration, cm^-3
};

// Per-material-pair segregation coefficients (P1-9), indexed by
// [min(matP,matN)][max(matP,matN)]. h is the interface transport rate
// (cm/s); m is the equilibrium ratio C_{lower id}/C_{higher id}. Defaults:
// h=0 (no exchange), m=1 (equal partition); the Si/oxide pair is filled by
// run()/run_ted() with segregation_h/m(dp, T) per P1-4.
struct SegTable {
  double h[5][5] = {};
  double m[5][5] = {{1, 1, 1, 1, 1}, {1, 1, 1, 1, 1}, {1, 1, 1, 1, 1},
                    {1, 1, 1, 1, 1}, {1, 1, 1, 1, 1}};
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
  // Per-cell material id (MatId: 0 Si, 1 oxide, 2 nitride, 3 poly, 4 gas).
  // The legacy value 2 (P1-4's "other/frozen") is synonymous with kMatGas.
  // Empty vector (default) derives ids from solve_mask: mask=1 -> kMatSi,
  // mask=0 -> kMatGas, i.e. the pre-P1-4 behavior.
  DiffusionSolver(const Mesh& mesh, std::vector<char> solve_mask,
                  std::ostream* log = nullptr,
                  std::vector<int> cell_mat = {});

  void run(std::vector<SpeciesField>& fields,
           const std::vector<DirichletBC>& bcs, const DiffuseOpts& opts);

  // Transient enhanced diffusion. Co-solves the excess self-interstitial field
  // `psi` (cm^-3 above equilibrium, seeded by the caller's "+1" implant model)
  // with the dopant fields. Each dopant's diffusivity is scaled by
  //   (1 - fi) + fi * (1 + psi / C_I*),
  // so the anneal starts with strongly enhanced diffusion that decays as the
  // interstitial excess diffuses to the surface sink and recombines.
  void run_ted(std::vector<SpeciesField>& fields, std::vector<double>& psi,
               const std::vector<DirichletBC>& bcs, const DiffuseOpts& opts);

 private:
  enum FaceKind : char {
    kInactive = 0,   // both sides masked
    kInternal = 1,   // flux between two active cells (same material)
    kBoundOwner = 2, // only owner active (boundary or material interface)
    kBoundNeigh = 3, // only neigh active (material interface)
    kSegregation = 4, // Si/oxide interface: segregation exchange
  };
  struct FGeom {
    char kind = kInactive;
    double g = 0;             // two-point geometric transmissibility |S|^2/(S.d)
    double wP = 0;            // face interpolation weight of the owner cell
    double delP = 0, delN = 0;  // centroid-to-face distances along the normal
    double gb = 0;            // boundary transmissibility |S|^2/(S.(cf-cP))
    Vec3 k;                   // non-orthogonal residual vector S - g*d
    double area = 0;          // |S|, kSegregation faces only
  };

  void build();
  void gradients(const std::vector<double>& c, const std::vector<double>& bcface,
                 std::vector<Vec3>& grad) const;

  // Assemble A_ and rhs for one backward-Euler scalar-transport step:
  //   (V/dt + reaction*V) c - div(D grad c) = (V/dt) cold   (+ nonortho corr)
  // `dcell` is the per-cell diffusivity, `bcface` the per-face Dirichlet value
  // (NaN = none), `cgrad` the field used for the deferred gradient correction
  // and to hold frozen (masked) cells. Uses the face coloring for parallelism.
  void assemble(const std::vector<double>& dcell,
                const std::vector<double>& cold,
                const std::vector<double>& bcface,
                const std::vector<double>& cgrad, double dt, double reaction,
                bool nonortho, const SegTable& seg,
                std::vector<double>& rhs, std::vector<Vec3>& grad);

  const Mesh& mesh_;
  std::vector<char> mask_;
  std::vector<int> mat_;   // per-cell MatId (0 Si .. 4 gas)
  bool has_segregation_ = false;  // any active kSegregation face
  // Distinct (lo,hi) MatId pairs present as kSegregation faces, filled by
  // build(); used to decide per-species whether SegTable has any m != 1
  // pair actually present in the mesh (-> switch to bicgstab_ilu0).
  std::vector<std::pair<int, int>> seg_pairs_;
  std::ostream* log_;
  std::vector<FGeom> fg_;
  CSR A_;
  std::vector<int> diag_;                  // diagonal slot per cell
  std::vector<std::array<int, 2>> fslot_;  // (P,N) and (N,P) slots per face
  int clamped_faces_ = 0;                  // faces with poor orthogonality

  // Greedy face coloring: faces within one color share no cell, so their
  // scatter-adds into cell rows are race-free and can run in parallel.
  std::vector<std::vector<int>> face_colors_;  // active face indices per color
};

}  // namespace cp
