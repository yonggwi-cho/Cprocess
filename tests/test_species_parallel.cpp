// PA-3: species-loop parallelization tests.
#include <cmath>
#include <cstdio>
#include <vector>

#include "cprocess/diffusion.hpp"
#include "cprocess/implant.hpp"
#include "cprocess/materials.hpp"
#include "cprocess/mesh.hpp"
#include "test_util.hpp"

using namespace cp;

namespace {

// 8x8x8 box mesh with B/P/As/Sb implanted near the top surface.
Mesh make_test_mesh() {
  const double um = 1e-4;
  return make_box_mesh(0, 1.0 * um, 0, 1.0 * um, 0, 1.0 * um, 8, 8, 8);
}

void implant_all(const Mesh& mesh, const std::vector<char>& mask,
                 std::vector<std::vector<double>>& conc,
                 const std::vector<const Dopant*>& dps) {
  const double um = 1e-4;
  conc.assign(dps.size(), std::vector<double>(mesh.cells.size(), 0.0));
  for (std::size_t k = 0; k < dps.size(); ++k) {
    ImplantParams ip;
    ip.dopant = dps[k];
    ip.dose = 1e14;
    ip.rp = 0.1 * um + 0.02 * um * k;
    ip.drp = 0.05 * um;
    apply_implant(mesh, mask, ip, conc[k]);
  }
}

}  // namespace

// 1. Equivalence: species_parallel = -1 (off) vs = 1 (on) must agree to
// within 1e-12 relative, cell by cell, for every species -- the per-species
// subproblems are fully independent so the parallel path uses the exact
// same math (assemble_into is byte-identical to assemble()'s former body).
static void test_equivalence() {
  Mesh mesh = make_test_mesh();
  std::vector<char> mask(mesh.cells.size(), 1);

  std::vector<const Dopant*> dps = {find_dopant("B"), find_dopant("P"),
                                    find_dopant("As"), find_dopant("Sb")};
  for (auto* d : dps) CHECK(d != nullptr);

  std::vector<std::vector<double>> conc_off, conc_on;
  implant_all(mesh, mask, conc_off, dps);
  conc_on = conc_off;

  DiffuseOpts o;
  o.temp = 1273.15;  // 1000 C
  o.time = 600;       // 10 min
  o.verbosity = 0;

  auto run_with = [&](std::vector<std::vector<double>>& conc, int sp) {
    DiffuseOpts oo = o;
    oo.species_parallel = sp;
    std::vector<SpeciesField> fields;
    for (std::size_t k = 0; k < dps.size(); ++k)
      fields.push_back({dps[k], &conc[k]});
    DiffusionSolver solver(mesh, mask, nullptr);
    solver.run(fields, {}, oo);
  };

  run_with(conc_off, -1);
  run_with(conc_on, 1);

  double maxrel = 0;
  for (std::size_t k = 0; k < dps.size(); ++k)
    for (std::size_t i = 0; i < conc_off[k].size(); ++i) {
      const double a = conc_off[k][i], b = conc_on[k][i];
      const double rel = std::fabs(a - b) / (std::fabs(a) + 1.0);
      maxrel = std::max(maxrel, rel);
    }
  std::printf("species_parallel equivalence: maxrel=%.3g\n", maxrel);
  CHECK(maxrel < 1e-12);
}

// 2. Auto heuristic with ns==1 must take the sp==false path (species_parallel
// left at its default 0): single-species diffuse should match the existing
// (non-species-parallel) diffusion tests bit-for-bit in structure -- here we
// simply check it runs and conserves dose, exercising the ns<3 auto branch.
static void test_auto_single_species() {
  Mesh mesh = make_test_mesh();
  std::vector<char> mask(mesh.cells.size(), 1);
  std::vector<const Dopant*> dps = {find_dopant("B")};
  std::vector<std::vector<double>> conc;
  implant_all(mesh, mask, conc, dps);

  double dose0 = 0;
  for (std::size_t i = 0; i < conc[0].size(); ++i)
    dose0 += conc[0][i] * mesh.cell_vol[i];

  DiffuseOpts o;  // species_parallel default 0 (auto); ns=1 -> sp==false
  o.temp = 1273.15;
  o.time = 300;
  o.verbosity = 0;
  std::vector<SpeciesField> fields = {{dps[0], &conc[0]}};
  DiffusionSolver solver(mesh, mask, nullptr);
  solver.run(fields, {}, o);

  double dose1 = 0;
  for (std::size_t i = 0; i < conc[0].size(); ++i) {
    CHECK(std::isfinite(conc[0][i]));
    dose1 += conc[0][i] * mesh.cell_vol[i];
  }
  CHECK_NEAR(dose1, dose0, 1e-3 * dose0);
}

// 3. Failure aggregation: force every species to fail to converge
// (lin_maxit=1) with species_parallel=1 -- the parallel path must not throw
// inside the OpenMP region (UB) but must still throw after it, once results
// are aggregated.
static void test_failure_aggregation() {
  Mesh mesh = make_test_mesh();
  std::vector<char> mask(mesh.cells.size(), 1);
  std::vector<const Dopant*> dps = {find_dopant("B"), find_dopant("P"),
                                    find_dopant("As")};
  std::vector<std::vector<double>> conc;
  implant_all(mesh, mask, conc, dps);

  DiffuseOpts o;
  o.temp = 1273.15;
  o.time = 60;
  o.dt = 60;
  o.verbosity = 0;
  o.species_parallel = 1;
  o.lin_maxit = 1;  // force non-convergence
  o.lin_rtol = 1e-14;

  std::vector<SpeciesField> fields;
  for (std::size_t k = 0; k < dps.size(); ++k) fields.push_back({dps[k], &conc[k]});
  DiffusionSolver solver(mesh, mask, nullptr);

  bool threw = false;
  try {
    solver.run(fields, {}, o);
  } catch (const std::exception&) {
    threw = true;
  }
  CHECK(threw);
}

int main() {
  test_equivalence();
  test_auto_single_species();
  test_failure_aggregation();
  std::printf("species_parallel tests passed\n");
  return 0;
}
