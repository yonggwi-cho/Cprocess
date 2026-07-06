#include <cmath>
#include <cstdio>
#include <vector>

#include "cprocess/diffusion.hpp"
#include "cprocess/implant.hpp"
#include "cprocess/materials.hpp"
#include "cprocess/mesh.hpp"
#include "test_util.hpp"

using namespace cp;

// Mass-weighted mean and variance of a profile along z.
static void moments(const Mesh& m, const std::vector<double>& c, double& mass,
                    double& mu, double& var) {
  mass = 0;
  double mz = 0;
  for (std::size_t i = 0; i < c.size(); ++i) {
    mass += c[i] * m.cell_vol[i];
    mz += c[i] * m.cell_vol[i] * m.cell_cent[i].z;
  }
  mu = mz / mass;
  double v = 0;
  for (std::size_t i = 0; i < c.size(); ++i)
    v += c[i] * m.cell_vol[i] * (m.cell_cent[i].z - mu) * (m.cell_cent[i].z - mu);
  var = v / mass;
}

// A buried Gaussian diffusing with constant D stays Gaussian with
// sigma^2(t) = sigma0^2 + 2 D t. Low concentration keeps D intrinsic.
static void test_gaussian_spread() {
  const double um = 1e-4;
  // Deliberately anisotropic cells (aspect ~5): exercises the deferred
  // non-orthogonal correction, without which the error here exceeds 30%.
  Mesh mesh = make_box_mesh(0, 0.2 * um, 0, 0.2 * um, 0, 2.0 * um, 2, 2, 96);

  const Dopant* boron = find_dopant("B");
  CHECK(boron != nullptr);

  ImplantParams ip;
  ip.dopant = boron;
  ip.dose = 1e12;        // peak ~8e16 cm^-3, well below ni(1373 K)
  ip.rp = 1.0 * um;      // buried: boundaries 7 sigma away after anneal
  ip.drp = 0.05 * um;
  std::vector<double> conc;
  std::vector<char> mask(mesh.cells.size(), 1);
  apply_implant(mesh, mask, ip, conc);

  double mass0, mu0, var0;
  moments(mesh, conc, mass0, mu0, var0);

  DiffuseOpts o;
  o.temp = 1373.15;  // 1100 C
  o.time = 600;
  o.dt = 30;
  o.field_enh = false;  // keep D strictly constant for the analytic check
  o.verbosity = 0;
  std::vector<SpeciesField> fields = {{boron, &conc}};
  DiffusionSolver solver(mesh, mask, nullptr);
  solver.run(fields, {}, o);

  double mass1, mu1, var1;
  moments(mesh, conc, mass1, mu1, var1);

  const double d = dopant_diffusivity(*boron, o.temp, 1.0);
  const double expect = 2.0 * d * o.time;
  std::printf("gaussian: D=%.4g cm^2/s  dvar=%.4g  expect=%.4g  (%.2f%%)\n", d,
              var1 - var0, expect, 100 * ((var1 - var0) / expect - 1));

  CHECK_NEAR(mass1, mass0, 1e-6 * mass0);          // zero-flux conservation
  CHECK_NEAR(var1 - var0, expect, 0.05 * expect);  // sigma^2 growth within 5%
  CHECK_NEAR(mu1, mu0, 0.01 * um);                 // no spurious drift
  for (double v : conc) CHECK(v >= 0.0);
}

// Constant-source predeposition: C(d,t) = Cs * erfc(d / (2 sqrt(D t))).
static void test_predeposition() {
  const double um = 1e-4;
  const double depth = 0.5 * um;
  Mesh mesh = make_box_mesh(0, 0.1 * um, 0, 0.1 * um, 0, depth, 2, 2, 100);

  const Dopant* boron = find_dopant("B");
  const double cs = 1e18;  // below ni(1223 K): near-intrinsic D

  DirichletBC bc;
  bc.species = "B";
  bc.patch = mesh.find_patch("zmax");
  bc.conc = cs;

  DiffuseOpts o;
  o.temp = 1223.15;  // 950 C
  o.time = 1200;
  o.dt = 60;
  o.field_enh = false;
  o.verbosity = 0;

  std::vector<double> conc(mesh.cells.size(), 0.0);
  std::vector<char> mask(mesh.cells.size(), 1);
  std::vector<SpeciesField> fields = {{boron, &conc}};
  DiffusionSolver solver(mesh, mask, nullptr);
  solver.run(fields, {bc}, o);

  const double d = dopant_diffusivity(*boron, o.temp, 1.0);
  const double l = 2.0 * std::sqrt(d * o.time);
  double err2 = 0, ref2 = 0;
  for (std::size_t i = 0; i < conc.size(); ++i) {
    const double dz = depth - mesh.cell_cent[i].z;
    const double ana = cs * std::erfc(dz / l);
    err2 += (conc[i] - ana) * (conc[i] - ana) * mesh.cell_vol[i];
    ref2 += ana * ana * mesh.cell_vol[i];
  }
  const double rel = std::sqrt(err2 / ref2);
  std::printf("predep: D=%.4g cm^2/s  2sqrt(Dt)=%.4g um  L2 rel err=%.3f\n", d,
              l / um, rel);
  CHECK(rel < 0.10);
  // Surface cells must approach Cs from below, never overshoot.
  for (double v : conc) CHECK(v <= 1.001 * cs);
}

// Masked (non-silicon) region: dopant must not leak through the interface.
static void test_mask_wall() {
  const double um = 1e-4;
  Mesh mesh = make_box_mesh(0, 0.2 * um, 0, 0.2 * um, 0, 1.0 * um, 2, 2, 20);
  std::vector<char> mask(mesh.cells.size(), 1);
  for (std::size_t i = 0; i < mesh.cells.size(); ++i)
    if (mesh.cell_cent[i].z > 0.5 * um) mask[i] = 0;  // top half frozen

  const Dopant* phos = find_dopant("P");
  std::vector<double> conc(mesh.cells.size(), 0.0);
  double active0 = 0;
  for (std::size_t i = 0; i < mesh.cells.size(); ++i)
    if (mask[i]) {
      conc[i] = 1e18;
      active0 += conc[i] * mesh.cell_vol[i];
    }

  DiffuseOpts o;
  o.temp = 1373.15;
  o.time = 600;
  o.dt = 60;
  o.verbosity = 0;
  std::vector<SpeciesField> fields = {{phos, &conc}};
  DiffusionSolver solver(mesh, mask, nullptr);
  solver.run(fields, {}, o);

  double active1 = 0;
  for (std::size_t i = 0; i < mesh.cells.size(); ++i) {
    if (mask[i]) active1 += conc[i] * mesh.cell_vol[i];
    else CHECK(conc[i] == 0.0);  // frozen cells untouched
  }
  CHECK_NEAR(active1, active0, 1e-6 * active0);
}

// PA-2: RCM reordering must not change physical results (dose conservation,
// reproducibility) or degrade convergence, relative to the pre-PA-2 solver.
static void test_rcm_integration() {
  const double um = 1e-4;
  Mesh mesh = make_box_mesh(0, 0.2 * um, 0, 0.2 * um, 0, 2.0 * um, 2, 2, 96);
  const Dopant* boron = find_dopant("B");
  CHECK(boron != nullptr);

  ImplantParams ip;
  ip.dopant = boron;
  ip.dose = 1e12;
  ip.rp = 1.0 * um;
  ip.drp = 0.05 * um;
  std::vector<double> conc0;
  std::vector<char> mask(mesh.cells.size(), 1);
  apply_implant(mesh, mask, ip, conc0);

  double mass0, mu0, var0;
  moments(mesh, conc0, mass0, mu0, var0);
  double peak0 = 0;
  for (double v : conc0) peak0 = std::max(peak0, v);

  DiffuseOpts o;
  o.temp = 1373.15;
  o.time = 600;
  o.dt = 30;
  o.field_enh = false;
  o.verbosity = 0;

  // --- run #1 ---
  std::vector<double> conc1 = conc0;
  int picard_iters1 = 0;
  o.nl_iters = &picard_iters1;
  std::vector<SpeciesField> fields1 = {{boron, &conc1}};
  DiffusionSolver solver1(mesh, mask, nullptr);
  solver1.run(fields1, {}, o);

  double mass1, mu1, var1;
  moments(mesh, conc1, mass1, mu1, var1);
  double peak1 = 0;
  for (double v : conc1) peak1 = std::max(peak1, v);

  // Measured on this exact case with the pre-PA-2 solver (commit e7568de,
  // natural mesh cell order, no RCM reordering): mass0=mass1=4.0e2 (exact,
  // as expected for zero-flux conservation), peak1=2.8170413934e16,
  // picard_iters=153. Baked in here as the reference for the "results
  // unchanged, convergence not degraded" checks below, per the spec's
  // prescribed method (measure once before the change, assert against the
  // recorded constant after).
  const double ref_mass = 4.0e2;
  const double ref_peak1 = 2.8170413934e16;
  const int ref_picard_iters = 153;

  std::printf("rcm integration: mass0=%.6e mass1=%.6e peak1=%.6e picard_iters=%d\n",
              mass0, mass1, peak1, picard_iters1);

  CHECK_NEAR(mass1, ref_mass, 1e-9 * ref_mass);  // dose conservation < 1e-9 rel
  CHECK_NEAR(peak1, ref_peak1, 1e-9 * ref_peak1);  // peak matches pre-PA-2
  CHECK(peak1 > 0 && peak1 < peak0);              // diffused, still positive

  // --- run #2: identical inputs must reproduce bit-identical results. ---
  std::vector<double> conc2 = conc0;
  int picard_iters2 = 0;
  o.nl_iters = &picard_iters2;
  std::vector<SpeciesField> fields2 = {{boron, &conc2}};
  DiffusionSolver solver2(mesh, mask, nullptr);
  solver2.run(fields2, {}, o);

  for (std::size_t i = 0; i < conc1.size(); ++i) CHECK(conc1[i] == conc2[i]);
  CHECK(picard_iters2 == picard_iters1);

  // --- iteration count within +/-10% of the pre-PA-2 baseline. ---
  CHECK(picard_iters1 <= static_cast<int>(1.10 * ref_picard_iters) + 1);
  CHECK(picard_iters1 >= static_cast<int>(0.90 * ref_picard_iters) - 1);
}

int main() {
  test_gaussian_spread();
  test_predeposition();
  test_mask_wall();
  test_rcm_integration();
  std::printf("diffusion tests passed\n");
  return 0;
}
