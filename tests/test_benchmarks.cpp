// Quantitative literature benchmarks — Tier A/C of the benchmark suite.
//
// Unlike test_sprocess_parity.cpp (directional parity vs the DESIRED
// Sentaurus-Process behavior, expected to fail until plan tasks land), this
// suite compares engine output against PUBLISHED, engine-independent
// reference data with generous tolerances, and must PASS:
//
//   Tier A — hard asserts. Engine agrees with the literature value within a
//            generous tolerance today; these are permanent regression
//            anchors tied to external truth.
//   Tier C — informational only ("INFO"). Either the reference value's
//            precision or the engine's model scope makes a hard assert
//            unjustified; measured-vs-reference is printed with the source
//            so drift is visible in the log, but nothing fails.
//
// Tier B (quantitative targets the engine currently misses by a large,
// unambiguous factor) lives in tests/test_sprocess_parity.cpp (WILL_FAIL),
// tagged with plan-task IDs.
//
// Data-integrity rule: every reference number below carries a source comment
// (author/year/textbook + conditions), and only very well-established,
// widely-reproduced values are used:
//  - Implant Rp/dRp: LSS-based range/straggle tables for B/P/As in Si as
//    reproduced in S.M. Sze, "Physics of Semiconductor Devices" (implant
//    appendix/Ch.14 tables) and Plummer/Deal/Griffin, "Silicon VLSI
//    Technology" Ch.8 (values rounded to the ~10% level those tables agree
//    to across editions).
//  - Oxidation: B.E. Deal & A.S. Grove, J. Appl. Phys. 36, 3770 (1965),
//    <100>-like rate constants as tabulated in the standard texts
//    (dry O2:  B=0.0117 um^2/h, B/A=0.071 um/h at 1000 C;
//              B=0.027  um^2/h, B/A=0.30  um/h at 1100 C, tau_dry~0.37 h
//              at 1000 C for the fast-initial-regime offset;
//     wet O2 (95 C H2O): B=0.287 um^2/h, B/A=1.27 um/h at 1000 C, tau~0).
//    Only the THICK regime (>50 nm) is asserted, where the Massoud
//    thin-oxide correction (absent from the engine, parity check C-3) does
//    not matter; thin-regime points are Tier C.
//  - Diffusivity: R.B. Fair's vacancy-model intrinsic boron diffusivity,
//    D_B = 0.76 exp(-3.46 eV/kT) cm^2/s ("Impurity Doping Processes in
//    Silicon", 1981; the classic SUPREM default), ~1.5e-14 cm^2/s at 1000 C.
//  - TED magnitude: classic marker-layer experiments (Packan & Plummer;
//    Stolk et al., J. Appl. Phys. 81, 6031 (1997)): time-averaged B
//    diffusivity enhancements of order 10-100x for short anneals at
//    750-810 C. Order-of-magnitude reference only -> Tier C here, hard band
//    in the parity suite (C-2).
//  - As electrical activation limit: ~2e20 cm^-3 at 900 C, ~3e20 cm^-3 at
//    1000 C (Nobili/Solmi electrical-activation studies, as summarized in
//    Plummer Ch.7). NOTE: the engine's own solid-solubility Arrhenius fit
//    derives from the same textbook data, so this comparison is partly
//    circular -> Tier C by design.

#include <cmath>
#include <cstdio>
#include <sstream>
#include <string>
#include <vector>

#include "cprocess/materials.hpp"
#include "cprocess/process.hpp"

using namespace cp;

// ---------------------------------------------------------------------------
// Harness
// ---------------------------------------------------------------------------

struct Row {
  std::string name;
  std::string source;
  double ref;        // reference value
  double meas;       // measured value
  double tol;        // relative tolerance (Tier A); informational for Tier C
  char tier;         // 'A' or 'C'
  bool pass;         // Tier A only; Tier C rows always "pass"
  std::string unit;
};

static std::vector<Row> g_rows;
static int g_fail = 0;

static void tier_a(const std::string& name, const std::string& source,
                   double ref, double meas, double tol,
                   const std::string& unit) {
  const bool ok = std::fabs(meas - ref) <= tol * std::fabs(ref);
  if (!ok) ++g_fail;
  g_rows.push_back({name, source, ref, meas, tol, 'A', ok, unit});
  std::printf("[A] %-34s ref=%-10.4g meas=%-10.4g tol=+-%.0f%%  %s\n",
              name.c_str(), ref, meas, tol * 100, ok ? "PASS" : "FAIL");
}

// Tier A asserted as a ratio band [lo, hi] (for factor-style tolerances).
static void tier_a_band(const std::string& name, const std::string& source,
                        double ref, double meas, double lo, double hi,
                        const std::string& unit) {
  const double r = meas / ref;
  const bool ok = r >= lo && r <= hi;
  if (!ok) ++g_fail;
  g_rows.push_back({name, source, ref, meas, hi - 1.0, 'A', ok, unit});
  std::printf("[A] %-34s ref=%-10.4g meas=%-10.4g band=[%.2g,%.2g]x  %s\n",
              name.c_str(), ref, meas, lo, hi, ok ? "PASS" : "FAIL");
}

static void tier_c(const std::string& name, const std::string& source,
                   double ref, double meas, const std::string& unit) {
  g_rows.push_back({name, source, ref, meas, 0, 'C', true, unit});
  std::printf("[C] %-34s ref=%-10.4g meas=%-10.4g  INFO\n", name.c_str(), ref,
              meas);
}

// ---------------------------------------------------------------------------
// Measurement helpers (fine 1-D column, proc:: calls only)
// ---------------------------------------------------------------------------

static SimState column(double zmax_cm, int nz) {
  SimState st;
  std::ostringstream l;
  proc::mesh_box(st, 0, 0.3e-4, 0, 0.3e-4, 0, zmax_cm, 2, 2, nz, &l);
  proc::set_region(st, "silicon", -1, &l);
  return st;
}

// Mass-weighted mean depth / spread and peak-cell depth of a species (cm).
struct Profile1D {
  double mean = 0, sigma = 0, peak = 0, dose_conc_sum = 0;
};

static Profile1D profile_moments(const SimState& st, const std::string& sp) {
  const auto& f = st.fields.at(sp);
  const double ztop = st.mesh.bbox().hi.z;
  double m = 0, md = 0, md2 = 0, best = -1;
  Profile1D r;
  for (std::size_t i = 0; i < f.size(); ++i) {
    const double w = f[i] * st.mesh.cell_vol[i];
    const double d = ztop - st.mesh.cell_cent[i].z;
    m += w;
    md += w * d;
    md2 += w * d * d;
    if (f[i] > best) {
      best = f[i];
      r.peak = d;
    }
  }
  r.dose_conc_sum = m;
  if (m > 0) {
    r.mean = md / m;
    r.sigma = std::sqrt(std::max(0.0, md2 / m - r.mean * r.mean));
  }
  return r;
}

// ---------------------------------------------------------------------------
// 1. Implant range/straggle: analytic Gaussian vs LSS/textbook tables
// ---------------------------------------------------------------------------
// Reference Rp/dRp (nm) in bare Si — LSS-based tables as reproduced in Sze,
// "Physics of Semiconductor Devices" and Plummer, "Silicon VLSI Technology"
// Ch.8 (values agree across editions to ~10%).
struct ImpRef {
  const char* sp;
  double e_kev, rp_nm, drp_nm;
};
static const ImpRef kImpRefs[] = {
    {"B", 30, 100, 34},   // B 30 keV:  Rp~0.100 um, dRp~0.034 um
    {"B", 100, 307, 69},  // B 100 keV: Rp~0.307 um, dRp~0.069 um
    {"P", 50, 61, 26},    // P 50 keV:  Rp~0.061 um, dRp~0.026 um
    {"P", 100, 123, 45},  // P 100 keV: Rp~0.123 um, dRp~0.045 um
    {"As", 50, 32, 12},   // As 50 keV: Rp~0.032 um, dRp~0.012 um
    {"As", 100, 58, 20},  // As 100 keV:Rp~0.058 um, dRp~0.020 um
};
static const char* kImpSrc = "LSS tables (Sze PSD; Plummer VLSI Ch.8)";

static void bench_implant_analytic() {
  for (const auto& r : kImpRefs) {
    SimState st = column(1.0e-4, 250);  // 4 nm cells
    std::ostringstream l;
    proc::implant_gauss(st, r.sp, 1e13, r.e_kev, 0, 0, 0, false, 0, 0, 0, 0,
                        false, "gauss", &l);
    const Profile1D p = profile_moments(st, r.sp);
    char n1[64], n2[64];
    std::snprintf(n1, sizeof n1, "analytic Rp %s %g keV", r.sp, r.e_kev);
    std::snprintf(n2, sizeof n2, "analytic dRp %s %g keV", r.sp, r.e_kev);
    // Measured today: B30 99.0/36.4, B100 299/71, P50 67/28.1, P100 125/44.4,
    // As50 34/12.8, As100 58/20.7 nm.
    tier_a(n1, kImpSrc, r.rp_nm, p.peak * 1e7, 0.25, "nm");
    tier_a(n2, kImpSrc, r.drp_nm, p.sigma * 1e7, 0.30, "nm");
  }
}

// ---------------------------------------------------------------------------
// 2. Implant range: Monte-Carlo BCA vs the SAME references. The MC transport
// is independent of the analytic moment table, so this genuinely
// cross-checks the BCA physics. channeling=false to match the amorphous-
// target LSS reference; fixed seed + threads=1 for determinism.
// Measured today (30k ions, mean depth): P50 71.9, P100 141.1, As50 38.8,
// As100 68.9 nm — within ~20% of LSS -> Tier A (+-30%). B is systematically
// ~35% deep (B30 133.5, B100 411.0 nm) — inside 1.5x, outside a defensible
// assert -> Tier C (light-ion electronic stopping is high in the BCA).
// ---------------------------------------------------------------------------
static void bench_implant_mc() {
  for (const auto& r : kImpRefs) {
    SimState st = column(1.0e-4, 250);
    std::ostringstream l;
    proc::implant_mc(st, r.sp, 1e13, r.e_kev, 30000, 0, 0, 7, 1,
                     /*channeling=*/false, false, 0, 0, 0, 0, false, false,
                     &l);
    const Profile1D p = profile_moments(st, r.sp);
    char n[64];
    std::snprintf(n, sizeof n, "MC Rp %s %g keV", r.sp, r.e_kev);
    if (std::string(r.sp) == "B")
      tier_c(n, kImpSrc, r.rp_nm, p.mean * 1e7, "nm");
    else
      tier_a(n, kImpSrc, r.rp_nm, p.mean * 1e7, 0.30, "nm");
  }
}

// ---------------------------------------------------------------------------
// 3. Oxidation vs Deal & Grove (1965). References computed from the D&G
// rate constants quoted in the header comment:
//   x(t) = A/2 * (sqrt(1 + (t+tau)/(A^2/4B)) - 1).
// Thick regime (>50 nm) -> Tier A; sub-60 nm dry points where the initial-
// regime offset tau (and the missing Massoud term, parity C-3) dominate the
// reference uncertainty -> Tier C.
// Measured today (proc::oxidize, default OED): dry1000C 30/60/120 min =
// 30.8/54.2/91.0 nm; dry1100C 30/60 = 74.8/117.7 nm; wet1000C 30/60 =
// 290.6/449.6 nm.
// ---------------------------------------------------------------------------
static const char* kDgSrc = "Deal & Grove, JAP 36, 3770 (1965)";

static double grow(double t_min, double temp_c, bool wet) {
  SimState st = column(0.8e-4, 40);
  std::ostringstream l;
  return proc::oxidize(st, t_min * 60.0, temp_c + 273.15, wet, &l) * 1e7;  // nm
}

static void bench_oxidation() {
  // Reference thicknesses (nm) from the D&G constants above:
  //  dry 1000 C, 120 min, tau=0.37 h : 103 nm
  //  dry 1100 C,  30 min, tau~0.10 h :  87 nm
  //  dry 1100 C,  60 min, tau~0.10 h : 130 nm
  //  wet 1000 C,  30 min             : 282 nm
  //  wet 1000 C,  60 min             : 435 nm
  tier_a("dry ox 1000C 120min", kDgSrc, 103, grow(120, 1000, false), 0.20,
         "nm");
  tier_a("dry ox 1100C 30min", kDgSrc, 87, grow(30, 1100, false), 0.25, "nm");
  tier_a("dry ox 1100C 60min", kDgSrc, 130, grow(60, 1100, false), 0.20, "nm");
  tier_a("wet ox 1000C 30min", kDgSrc, 282, grow(30, 1000, true), 0.15, "nm");
  tier_a("wet ox 1000C 60min", kDgSrc, 435, grow(60, 1000, true), 0.15, "nm");
  // Thin-regime dry points: reference itself is tau/Massoud-limited.
  tier_c("dry ox 1000C 30min (thin)", kDgSrc, 48, grow(30, 1000, false), "nm");
  tier_c("dry ox 1000C 60min (thin)", kDgSrc, 69, grow(60, 1000, false), "nm");
}

// ---------------------------------------------------------------------------
// 4. Boron intrinsic diffusivity: buried Gaussian marker (peak 2e18 cm^-3 <
// n_i(1273 K) ~ 7e18, i.e. intrinsic), 1000 C / 1 h drive. The variance of a
// Gaussian grows as sigma^2(t) = sigma0^2 + 2 D t, so
// (sigma1^2 - sigma0^2)/(2 t) is a direct diffusivity measurement compared
// against Fair's D_B = 0.76 exp(-3.46 eV/kT) cm^2/s. Buried at 1 um so the
// surface boundary never truncates the spread. Factor-2 band (task-level
// tolerance for diffusivity-driven motion); measured today: ratio 1.03.
// ---------------------------------------------------------------------------
static const char* kFairSrc = "Fair 1981 (D_B=0.76*exp(-3.46eV/kT))";

static void bench_boron_diffusivity() {
  SimState st = column(2.0e-4, 400);  // 5 nm cells
  std::ostringstream l;
  proc::implant_gauss(st, "B", 1e13, 0, 1.0e-4, 0.02e-4, 0, false, 0, 0, 0, 0,
                      false, "gauss", &l);
  const double s0 = profile_moments(st, "B").sigma;
  DiffuseOpts d;
  d.temp = 1273.15;
  d.time = 3600;
  d.verbosity = 0;
  proc::diffuse(st, d, &l);
  const double s1 = profile_moments(st, "B").sigma;
  const double d_meas = (s1 * s1 - s0 * s0) / (2.0 * d.time);
  const double kt = 8.617333e-5 * d.temp;
  const double d_lit = 0.76 * std::exp(-3.46 / kt);  // ~1.5e-14 cm^2/s
  tier_a_band("B intrinsic D, 1000C marker", kFairSrc, d_lit, d_meas, 0.5,
              2.0, "cm2/s");
}

// ---------------------------------------------------------------------------
// 5. Junction depth of an implant+drive vs the closed-form Gaussian solution
// with Fair's literature D. B 1e14 cm^-2 (Rp 50 nm, sigma0 20 nm), drive
// 1100 C / 30 min, junction against a 1e15 cm^-3 background:
//   sigma^2 = sigma0^2 + 2 D t,  Cpk = Q/sqrt(2 pi sigma^2),
//   xj = Rp + sigma * sqrt(2 ln(Cpk/CB)).
// Peak stays < n_i(1373 K) ~ 1.2e19 -> intrinsic. Measured today:
// engine xj = 1.008 um vs analytic 0.954 um (+5.7%).
// ---------------------------------------------------------------------------
static void bench_junction_depth() {
  SimState st = column(2.0e-4, 200);
  std::ostringstream l;
  proc::implant_gauss(st, "B", 1e14, 0, 0.05e-4, 0.02e-4, 0, false, 0, 0, 0, 0,
                      false, "gauss", &l);
  DiffuseOpts d;
  d.temp = 1373.15;
  d.time = 1800;
  d.verbosity = 0;
  proc::diffuse(st, d, &l);
  const auto& f = st.fields.at("B");
  const double ztop = st.mesh.bbox().hi.z;
  double xj = 0;
  for (std::size_t i = 0; i < f.size(); ++i)
    if (f[i] >= 1e15) xj = std::max(xj, ztop - st.mesh.cell_cent[i].z);
  const double kt = 8.617333e-5 * d.temp;
  const double d_lit = 0.76 * std::exp(-3.46 / kt);
  const double var = 0.02e-4 * 0.02e-4 + 2 * d_lit * d.time;
  const double cpk = 1e14 / std::sqrt(2 * M_PI * var);
  const double xja = 0.05e-4 + std::sqrt(2 * var * std::log(cpk / 1e15));
  tier_a("B drive-in xj (1100C 30min)", kFairSrc, xja * 1e4, xj * 1e4, 0.30,
         "um");
}

// ---------------------------------------------------------------------------
// 6. TED enhancement magnitude (Tier C here; hard band lives in the parity
// suite, C-2). B marker + damage, 900 C / 60 s; time-averaged Dt enhancement
// (sigma_ted^2 - sigma0^2)/(sigma_eq^2 - sigma0^2) vs the classic 10-100x
// marker-experiment range (Packan/Plummer; Stolk 1997 — measured at
// 750-810 C, order-of-magnitude reference only). Measured today: ~706x.
// (The 800 C case currently produces NaN — parity check C-2.)
// ---------------------------------------------------------------------------
static void bench_ted() {
  DiffuseOpts d;
  d.temp = 1173.15;
  d.time = 60;
  d.verbosity = 0;
  std::ostringstream l;
  SimState eq = column(1.0e-4, 100);
  proc::implant_gauss(eq, "B", 1e14, 0, 0.05e-4, 0.02e-4, 0, false, 0, 0, 0, 0,
                      false, "gauss", &l);
  const double s0 = profile_moments(eq, "B").sigma;
  proc::diffuse(eq, d, &l);
  const double seq = profile_moments(eq, "B").sigma;
  SimState td = column(1.0e-4, 100);
  proc::implant_gauss(td, "B", 1e14, 0, 0.05e-4, 0.02e-4, 0, false, 0, 0, 0, 0,
                      true, "gauss", &l);
  proc::diffuse_ted(td, d, &l);
  const double sted = profile_moments(td, "B").sigma;
  const double enh = (sted * sted - s0 * s0) / (seq * seq - s0 * s0);
  tier_c("TED Dt-enhancement 900C 60s",
         "Packan/Stolk marker expts: 10-100x @750-810C", 30, enh, "x");
}

// ---------------------------------------------------------------------------
// 7. As electrical-activation clamp vs literature (~2e20 @900 C, ~3e20
// @1000 C; Nobili/Solmi as summarized in Plummer Ch.7). Tier C: the engine's
// own ss Arrhenius fit derives from the same textbook data, so a hard assert
// would be partly circular. Measured today: 1.90e20 / 3.17e20 cm^-3.
// ---------------------------------------------------------------------------
static void bench_as_activation() {
  SimState st = column(0.5e-4, 25);
  std::ostringstream l;
  proc::init(st, "As", 1e21, -1, &l);
  auto clamp = [&](double temp_c) {
    auto a = proc::active_field(st, "As", temp_c + 273.15, &l);
    double m = 0;
    for (double v : a) m = std::max(m, v);
    return m;
  };
  const char* src = "Nobili/Solmi (Plummer Ch.7); circular w/ engine fit";
  tier_c("As active limit 900C", src, 2e20, clamp(900), "cm-3");
  tier_c("As active limit 1000C", src, 3e20, clamp(1000), "cm-3");
}

// ---------------------------------------------------------------------------
int main() {
  std::printf("Quantitative literature benchmarks (Tier A assert / Tier C info)\n");
  std::printf("----------------------------------------------------------------\n");
  bench_implant_analytic();
  bench_implant_mc();
  bench_oxidation();
  bench_boron_diffusivity();
  bench_junction_depth();
  bench_ted();
  bench_as_activation();

  std::printf("\n============================ benchmark summary ============================\n");
  std::printf("%-4s %-34s %-12s %-12s %-8s %s\n", "tier", "check", "reference",
              "measured", "status", "source");
  int na = 0, nc = 0;
  for (const auto& r : g_rows) {
    std::printf("%-4c %-34s %-12.4g %-12.4g %-8s %s\n", r.tier, r.name.c_str(),
                r.ref, r.meas, r.tier == 'C' ? "INFO" : (r.pass ? "PASS" : "FAIL"),
                r.source.c_str());
    if (r.tier == 'A') ++na;
    else ++nc;
  }
  std::printf("===========================================================================\n");
  std::printf("%d Tier A checks (%d failed), %d Tier C informational rows\n", na,
              g_fail, nc);
  if (g_fail) {
    std::printf("BENCHMARK FAILURES: %d\n", g_fail);
    return 1;
  }
  std::printf("all Tier A benchmarks pass\n");
  return 0;
}
