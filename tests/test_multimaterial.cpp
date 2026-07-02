#include <cmath>
#include <cstdio>
#include <vector>

#include "cprocess/diffusion.hpp"
#include "cprocess/materials.hpp"
#include "cprocess/mesh.hpp"
#include "test_util.hpp"

using namespace cp;

namespace {

const double kUm = 1e-4;

// Box mesh: 0.05x0.05x0.4 um, 2x2x40 cells (cell height 0.01 um).
Mesh make_test_mesh() { return make_box_mesh(0, 0.05 * kUm, 0, 0.05 * kUm, 0, 0.4 * kUm, 2, 2, 40); }

double sum_cv(const Mesh& m, const std::vector<double>& c) {
  double s = 0;
  for (std::size_t i = 0; i < m.cells.size(); ++i) s += c[i] * m.cell_vol[i];
  return s;
}

// Gaussian implant-like profile in z, centered at z0 (cm) with sigma (cm),
// peak conc, placed only in cells matching `predicate`.
template <typename Pred>
std::vector<double> gaussian_field(const Mesh& m, double z0, double sigma,
                                   double peak, Pred pred) {
  std::vector<double> c(m.cells.size(), 0.0);
  for (std::size_t i = 0; i < m.cells.size(); ++i) {
    if (!pred(i)) continue;
    const double z = m.cell_cent[i].z;
    c[i] = peak * std::exp(-0.5 * (z - z0) * (z - z0) / (sigma * sigma));
  }
  return c;
}

// 1. Nitride barrier: top 0.1 um nitride, rest Si. B Gaussian in Si only.
void test_nitride_barrier() {
  Mesh m = make_test_mesh();
  std::vector<int> mat(m.cells.size(), kMatSi);
  for (std::size_t i = 0; i < m.cells.size(); ++i)
    if (m.cell_cent[i].z > 0.3 * kUm) mat[i] = kMatNitride;

  const Dopant* boron = find_dopant("B");
  CHECK(boron != nullptr);
  // Rp ~ 0.2 um deep in the Si region, comfortably below the interface.
  std::vector<double> c = gaussian_field(m, 0.2 * kUm, 0.03 * kUm, 1e19,
                                         [&](std::size_t i) { return mat[i] == kMatSi; });
  const double total0 = sum_cv(m, c);

  DiffuseOpts o;
  o.temp = 1273.15;  // 1000 C
  o.time = 1800;      // 30 min
  o.dt = 30;
  o.verbosity = 0;
  std::vector<SpeciesField> fields = {{boron, &c}};
  DiffusionSolver solver(m, {}, nullptr, mat);
  solver.run(fields, {}, o);

  double nit_dose = 0, total1 = 0;
  for (std::size_t i = 0; i < m.cells.size(); ++i) {
    total1 += c[i] * m.cell_vol[i];
    if (mat[i] == kMatNitride) nit_dose += c[i] * m.cell_vol[i];
  }
  std::printf("nitride barrier: nit_dose/total=%.4g%%  consv=%.4g%%\n",
              100 * nit_dose / total1, 100 * (total1 - total0) / total0);
  CHECK(nit_dose / total1 < 0.001);
  CHECK_NEAR(total1, total0, 0.005 * total0);
}

// 2. Poly fast path: spread-growth ratio vs. all-Si column.
void test_poly_fast_path() {
  const double z0 = 0.2 * kUm, sigma = 0.02 * kUm, peak = 1e18;

  auto run_column = [&](int mat_id, double& var0, double& var1) {
    Mesh m = make_test_mesh();
    std::vector<int> mat(m.cells.size(), mat_id);
    const Dopant* boron = find_dopant("B");
    std::vector<double> c = gaussian_field(m, z0, sigma, peak, [](std::size_t) { return true; });

    auto variance = [&](const std::vector<double>& c) {
      double w = 0, mean = 0;
      for (std::size_t i = 0; i < m.cells.size(); ++i) {
        w += c[i] * m.cell_vol[i];
        mean += c[i] * m.cell_vol[i] * m.cell_cent[i].z;
      }
      mean /= w;
      double v = 0;
      for (std::size_t i = 0; i < m.cells.size(); ++i)
        v += c[i] * m.cell_vol[i] * (m.cell_cent[i].z - mean) * (m.cell_cent[i].z - mean);
      return v / w;
    };
    var0 = variance(c);

    // A shorter anneal than the other tests: with 30 min at 1000 C the poly
    // column's spread (D_poly ~ 10x Si) becomes comparable to the distance
    // to the zero-flux top/bottom walls, and reflection compresses the
    // measured variance growth well below the free-diffusion ratio. 5 min
    // keeps both columns' spread << domain half-height (0.2 um) so the
    // ratio reflects the diffusivity ratio, not wall confinement.
    DiffuseOpts o;
    o.temp = 1273.15;
    o.time = 300;
    o.dt = 6;
    o.field_enh = false;
    o.verbosity = 0;
    std::vector<SpeciesField> fields = {{boron, &c}};
    DiffusionSolver solver(m, {}, nullptr, mat);
    solver.run(fields, {}, o);
    var1 = variance(c);
  };

  double var0_si, var1_si, var0_poly, var1_poly;
  run_column(kMatSi, var0_si, var1_si);
  run_column(kMatPoly, var0_poly, var1_poly);

  const double ratio = (var1_poly - var0_poly) / (var1_si - var0_si);
  std::printf("poly fast path: dvar_si=%.4g dvar_poly=%.4g ratio=%.4g\n",
              var1_si - var0_si, var1_poly - var0_poly, ratio);
  CHECK(ratio >= 5.0 && ratio <= 20.0);
}

// 3. Si-only regression: cell_mat all-Si vs. mask-derived solver agree tightly.
void test_si_only_regression() {
  Mesh m = make_test_mesh();
  const Dopant* boron = find_dopant("B");
  const Dopant* phos = find_dopant("P");

  std::vector<int> mat_all_si(m.cells.size(), kMatSi);
  std::vector<char> mask_all(m.cells.size(), 1);

  std::vector<double> b_new(m.cells.size(), 1e18), p_new(m.cells.size(), 5e17);
  std::vector<double> b_old(m.cells.size(), 1e18), p_old(m.cells.size(), 5e17);

  DiffuseOpts o;
  o.temp = 1273.15;
  o.time = 600;
  o.dt = 30;
  o.field_enh = true;
  o.verbosity = 0;

  std::vector<SpeciesField> fields_new = {{boron, &b_new}, {phos, &p_new}};
  DiffusionSolver solver_new(m, {}, nullptr, mat_all_si);
  solver_new.run(fields_new, {}, o);

  std::vector<SpeciesField> fields_old = {{boron, &b_old}, {phos, &p_old}};
  DiffusionSolver solver_old(m, mask_all, nullptr);
  solver_old.run(fields_old, {}, o);

  double maxrel = 0;
  double dose_new = sum_cv(m, b_new) + sum_cv(m, p_new);
  double dose_old = sum_cv(m, b_old) + sum_cv(m, p_old);
  for (std::size_t i = 0; i < m.cells.size(); ++i) {
    maxrel = std::max(maxrel, std::fabs(b_new[i] - b_old[i]) / (std::fabs(b_old[i]) + 1.0));
    maxrel = std::max(maxrel, std::fabs(p_new[i] - p_old[i]) / (std::fabs(p_old[i]) + 1.0));
  }
  std::printf("si-only regression: maxrel=%.4g dose_reldiff=%.4g\n", maxrel,
              std::fabs(dose_new - dose_old) / dose_old);
  CHECK(maxrel < 1e-6);
}

// 4. Multi-material stack dose conservation + m=1 equal-partition at Si/poly.
void test_stack_conservation() {
  Mesh m = make_test_mesh();
  std::vector<int> mat(m.cells.size(), kMatSi);
  for (std::size_t i = 0; i < m.cells.size(); ++i) {
    const double z = m.cell_cent[i].z;
    if (z > 0.3 * kUm) mat[i] = kMatOxide;        // top 0.1 um
    else if (z > 0.2 * kUm) mat[i] = kMatPoly;    // middle 0.1 um
    // bottom 0.2 um stays Si
  }
  const Dopant* boron = find_dopant("B");
  std::vector<double> c(m.cells.size(), 0.0);
  for (std::size_t i = 0; i < m.cells.size(); ++i)
    if (mat[i] == kMatSi || mat[i] == kMatPoly) c[i] = 1e18;
  const double total0 = sum_cv(m, c);

  DiffuseOpts o;
  o.temp = 1273.15;
  o.time = 1800;
  o.dt = 30;
  o.verbosity = 0;
  std::vector<SpeciesField> fields = {{boron, &c}};
  DiffusionSolver solver(m, {}, nullptr, mat);
  solver.run(fields, {}, o);

  const double total1 = sum_cv(m, c);
  std::printf("stack conservation: total0=%.6g total1=%.6g reldiff=%.4g%%\n",
              total0, total1, 100 * (total1 - total0) / total0);
  CHECK_NEAR(total1, total0, 0.005 * total0);

  // Adjacent Si/poly interface cells: equal partition (m=1) within 10%.
  double si_z = -1e30, poly_z = 1e30;
  for (std::size_t i = 0; i < m.cells.size(); ++i) {
    if (mat[i] == kMatSi) si_z = std::max(si_z, m.cell_cent[i].z);
    if (mat[i] == kMatPoly) poly_z = std::min(poly_z, m.cell_cent[i].z);
  }
  double si_sum = 0, poly_sum = 0;
  int si_n = 0, poly_n = 0;
  for (std::size_t i = 0; i < m.cells.size(); ++i) {
    if (mat[i] == kMatSi && std::fabs(m.cell_cent[i].z - si_z) < 1e-9) { si_sum += c[i]; ++si_n; }
    if (mat[i] == kMatPoly && std::fabs(m.cell_cent[i].z - poly_z) < 1e-9) { poly_sum += c[i]; ++poly_n; }
  }
  const double ratio = (poly_sum / poly_n) / (si_sum / si_n);
  std::printf("stack Si/poly interface ratio=%.4g\n", ratio);
  CHECK(ratio > 0.9 && ratio < 1.1);
}

}  // namespace

int main() {
  test_nitride_barrier();
  test_poly_fast_path();
  test_si_only_regression();
  test_stack_conservation();
  std::printf("multimaterial tests passed\n");
  return 0;
}
