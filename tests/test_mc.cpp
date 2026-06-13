#include <cmath>
#include <cstdio>
#include <vector>

#include "cprocess/cell_locator.hpp"
#include "cprocess/mc_implant.hpp"
#include "cprocess/mesh.hpp"
#include "test_util.hpp"

using namespace cp;

// ---------------------------------------------------------------------------
// Reference: exact classical scattering integral for the ZBL potential.
//   theta = pi - 2 I,  I = int_{R0}^inf  b dR / (R^2 sqrt(F(R))),
//   F(R) = 1 - V(R)/eps - b^2/R^2,  V = Phi(R)/R.
// Substituting R = R0/(1 - w^2) removes the sqrt singularity at R0.
// Then sin^2(theta/2) = cos^2(I).
// ---------------------------------------------------------------------------
static double zbl_phi(double x) {
  static const double c[4] = {0.18175, 0.50986, 0.28022, 0.02817};
  static const double d[4] = {-3.19980, -0.94229, -0.40290, -0.20162};
  double s = 0;
  for (int i = 0; i < 4; ++i) s += c[i] * std::exp(d[i] * x);
  return s;
}

static double exact_sin2half(double eps, double b) {
  auto F = [&](double r) {
    return 1.0 - zbl_phi(r) / (eps * r) - b * b / (r * r);
  };
  // F is monotonically increasing; bracket and bisect for R0.
  double hi = std::max(b, 1.0);
  while (F(hi) < 0) hi *= 2;
  double lo = hi * 0.5;
  while (lo > 1e-14 && F(lo) > 0) { hi = lo; lo *= 0.5; }
  for (int i = 0; i < 200; ++i) {
    const double mid = 0.5 * (lo + hi);
    (F(mid) > 0 ? hi : lo) = mid;
  }
  const double r0 = 0.5 * (lo + hi);

  // dF/dR at R0 for the analytic limit of the integrand at w -> 0.
  const double h = 1e-6 * r0;
  const double dF = (F(r0 + h) - F(r0 - h)) / (2 * h);
  const double lim0 = 2.0 * b / (r0 * std::sqrt(std::max(dF * r0, 1e-300)));

  // With R = r0/u and u = 1-w^2:
  //   I = int b/(R^2 sqrt(F)) dR = int_0^1 (b/r0) du / sqrt(F(r0/u))
  //     = int_0^1 (b/r0) 2 w dw / sqrt(F(r0/(1-w^2))),
  // integrated with composite Simpson; near w=0 the analytic limit avoids
  // the 0/0 from catastrophic cancellation in F(R0).
  const int n = 4000;  // even
  auto g = [&](double wv) {
    if (wv < 1e-4) return lim0;
    const double r = r0 / (1.0 - wv * wv);
    return 2.0 * wv * (b / r0) / std::sqrt(std::max(F(r), 1e-300));
  };
  double sum = g(0.0) + g(1.0);
  for (int i = 1; i < n; ++i) sum += g(static_cast<double>(i) / n) * (i % 2 ? 4.0 : 2.0);
  const double I = sum / (3.0 * n);

  const double c = std::cos(I);
  return c * c;
}

static void test_table_vs_integral() {
  const double epss[] = {0.001, 0.01, 0.1, 1.0, 10.0, 100.0};
  const double bs[] = {0.05, 0.3, 1.0, 2.0, 5.0};
  for (double eps : epss)
    for (double b : bs) {
      const double ex = exact_sin2half(eps, b);
      const double tb = mc::sin2_half_theta(eps, b);
      const double tol = std::max(0.02 * ex, 2e-4);
      std::printf("scatter eps=%-6g b=%-4g exact=%.6f table=%.6f\n", eps, b,
                  ex, tb);
      CHECK_NEAR(tb, ex, tol);
    }
}

// ---------------------------------------------------------------------------
static McImplantStats run_implant(const Mesh& mesh, long long ions,
                                  int threads, std::uint64_t seed,
                                  std::vector<double>& conc) {
  McImplantParams p;
  p.dopant = find_dopant("B");
  p.dose = 1e13;
  p.energy_kev = 50;
  p.ions = ions;
  p.threads = threads;
  p.seed = seed;
  std::vector<char> mask(mesh.cells.size(), 1);
  conc.assign(mesh.cells.size(), 0.0);
  return apply_mc_implant(mesh, mask, p, conc);
}

static void test_physics_and_binning() {
  const double um = 1e-4;
  Mesh mesh = make_box_mesh(0, 0.4 * um, 0, 0.4 * um, 0, 0.8 * um, 8, 8, 80);

  std::vector<double> conc;
  const long long ions = 20000;
  const McImplantStats s = run_implant(mesh, ions, 0, 42, conc);
  std::printf("B 50 keV: Rp=%.1f nm dRp=%.1f nm dep=%lld back=%lld\n",
              s.rp * 1e7, s.drp * 1e7, s.deposited, s.backscattered);

  // Accounting closes exactly.
  CHECK(s.deposited + s.backscattered + s.transmitted + s.out_of_domain +
            s.in_mask + s.unbinned == ions);
  CHECK(s.unbinned == 0);
  CHECK(s.transmitted == 0);          // 0.8 um box swallows 50 keV boron
  CHECK(s.out_of_domain == 0);        // full-area implant wraps laterally
  CHECK(s.backscattered < ions / 20); // a few % at most for vertical B

  // Physical ballpark for B 50 keV in amorphous Si (tables: ~161 nm / 50 nm).
  CHECK(s.rp > 100e-7 && s.rp < 250e-7);
  CHECK(s.drp > 30e-7 && s.drp < 95e-7);

  // The binned field integrates to exactly weight * deposited.
  double field_atoms = 0, mz = 0;
  for (std::size_t i = 0; i < conc.size(); ++i) {
    field_atoms += conc[i] * mesh.cell_vol[i];
    mz += conc[i] * mesh.cell_vol[i] * (mesh.bbox().hi.z - mesh.cell_cent[i].z);
  }
  CHECK_NEAR(field_atoms, s.weight * s.deposited, 1e-9 * field_atoms);
  // Binning preserves the depth moment (cell-centroid quantization only).
  CHECK_NEAR(mz / field_atoms, s.rp, 0.03 * s.rp);
}

static void test_thread_determinism() {
  const double um = 1e-4;
  Mesh mesh = make_box_mesh(0, 0.3 * um, 0, 0.3 * um, 0, 0.6 * um, 6, 6, 40);
  std::vector<double> c1, c4;
  const McImplantStats s1 = run_implant(mesh, 6000, 1, 7, c1);
  const McImplantStats s4 = run_implant(mesh, 6000, 4, 7, c4);
  CHECK(s1.deposited == s4.deposited);
  CHECK(s1.backscattered == s4.backscattered);
  CHECK(s1.rp == s4.rp);  // bit-identical reduction
  for (std::size_t i = 0; i < c1.size(); ++i) CHECK(c1[i] == c4[i]);
  std::printf("determinism: 1 vs 4 threads bit-identical (dep=%lld)\n",
              s1.deposited);
}

static void test_window() {
  const double um = 1e-4;
  Mesh mesh = make_box_mesh(0, 0.6 * um, 0, 0.6 * um, 0, 0.3 * um, 12, 12, 24);
  McImplantParams p;
  p.dopant = find_dopant("As");
  p.dose = 1e14;
  p.energy_kev = 30;
  p.ions = 6000;
  p.seed = 3;
  p.has_window = true;
  p.x1 = 0.2 * um; p.x2 = 0.4 * um;
  p.y1 = 0.2 * um; p.y2 = 0.4 * um;
  std::vector<char> mask(mesh.cells.size(), 1);
  std::vector<double> conc(mesh.cells.size(), 0.0);
  const McImplantStats s = apply_mc_implant(mesh, mask, p, conc);

  CHECK(s.deposited > 0.9 * p.ions);  // shallow As: nearly all retained
  // Lateral straggle of 30 keV As is ~10 nm; far outside the window the
  // concentration must vanish.
  double far = 0, peak = 0;
  for (std::size_t i = 0; i < conc.size(); ++i) {
    peak = std::max(peak, conc[i]);
    const Vec3& cc = mesh.cell_cent[i];
    if (cc.x < 0.1 * um || cc.x > 0.5 * um || cc.y < 0.1 * um ||
        cc.y > 0.5 * um)
      far = std::max(far, conc[i]);
  }
  CHECK(peak > 0);
  CHECK(far < 0.01 * peak);
  std::printf("window: peak=%.3g cm^-3, far-field/peak=%.2g\n", peak,
              far / std::max(peak, 1.0));
}

// ---------------------------------------------------------------------------
// Channeling test: B 50 keV along Si <100> should have deeper Rp than the
// amorphous case (channeling "tail" shifts the mean depth).  We also verify
// that enabling channeling does not break dose accounting.
// ---------------------------------------------------------------------------
static void test_channeling() {
  const double um = 1e-4;
  // Deeper box (1.5 um) so channeled ions don't transmit.
  Mesh mesh = make_box_mesh(0, 0.4*um, 0, 0.4*um, 0, 1.5*um, 8, 8, 120);

  auto run = [&](bool ch) {
    McImplantParams p;
    p.dopant = find_dopant("B");
    p.dose   = 1e13;
    p.energy_kev = 50;
    p.ions   = 20000;
    p.seed   = 17;
    p.channeling = ch;
    std::vector<char> mask(mesh.cells.size(), 1);
    std::vector<double> conc(mesh.cells.size(), 0.0);
    return apply_mc_implant(mesh, mask, p, conc);
  };

  const McImplantStats sa = run(false);
  const McImplantStats sc = run(true);

  std::printf("channeling: amorphous Rp=%.1f nm, crystal Rp=%.1f nm\n",
              sa.rp*1e7, sc.rp*1e7);

  // Channeling ions travel deeper on average (Rp_crystal > Rp_amorphous).
  CHECK(sc.rp > sa.rp);
  // Dose accounting must still close.
  const long long ions = 20000;
  CHECK(sc.deposited + sc.backscattered + sc.transmitted +
        sc.out_of_domain + sc.in_mask + sc.unbinned == ions);
}

int main() {
  test_table_vs_integral();
  test_physics_and_binning();
  test_thread_determinism();
  test_window();
  test_channeling();
  std::printf("mc tests passed\n");
  return 0;
}
