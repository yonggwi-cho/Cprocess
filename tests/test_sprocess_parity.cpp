// SProcess-parity executable specification (Tier 2 of the
// IMPLEMENTATION_PLAN_v2 test specification).
//
// Every check below asserts the DESIRED (Sentaurus-Process-correct) behavior
// using only currently-existing proc:: APIs, so this file always compiles —
// but the checks are EXPECTED TO FAIL until the corresponding plan task
// lands. The CTest registration therefore uses WILL_FAIL TRUE: the suite
// stays green while gaps remain. When a fix lands and a check starts
// passing, ctest reports this test as FAILED — that is the signal to move
// the now-passing check into the regular suites (test_golden_flows /
// test_<feature>) and keep the remainder here.
//
// The harness runs ALL checks (no early exit), prints a per-check PASS/FAIL
// summary table tagged with the IMPLEMENTATION_PLAN_v2 task ID plus the
// measured values, and exits nonzero if ANY check fails.
//
// Current measured behavior (2026-07, HEAD 86d768f) is quoted per check so a
// future regression in the *measurement setup* (as opposed to the physics)
// is recognizable.

#include <cmath>
#include <cstdio>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

#include "cprocess/diffusion.hpp"
#include "cprocess/oxidation.hpp"
#include "cprocess/process.hpp"

using namespace cp;

// ---------------------------------------------------------------------------
// Harness
// ---------------------------------------------------------------------------

struct CheckResult {
  std::string task;   // plan task ID, e.g. "W-8"
  std::string name;
  bool pass = false;
  std::string detail; // measured values
};

static std::vector<CheckResult> g_results;

static void record(const std::string& task, const std::string& name,
                   bool pass, const std::string& detail) {
  g_results.push_back({task, name, pass, detail});
  std::printf("[%-6s] %-38s %s  (%s)\n", task.c_str(), name.c_str(),
              pass ? "PASS" : "FAIL", detail.c_str());
}

static std::string fmtv(const char* f, double a, double b) {
  char buf[160];
  std::snprintf(buf, sizeof buf, f, a, b);
  return buf;
}

static std::string mat_of(const SimState& st, int ci) {
  auto it = st.region_material.find(st.mesh.cell_region[ci]);
  return it == st.region_material.end() ? std::string("silicon") : it->second;
}

// Peak-concentration depth of `sp` below z_ref, silicon cells only (cm).
static double si_peak_depth(const SimState& st, const std::string& sp,
                            double z_ref) {
  const auto& f = st.fields.at(sp);
  double best = -1, zbest = z_ref;
  for (std::size_t i = 0; i < f.size(); ++i) {
    if (mat_of(st, (int)i) != "silicon") continue;
    if (f[i] > best) { best = f[i]; zbest = st.mesh.cell_cent[i].z; }
  }
  return z_ref - zbest;
}

static SimState make_column(double zmax, int nz, std::ostream* log) {
  SimState st;
  proc::mesh_box(st, 0, 0.3e-4, 0, 0.3e-4, 0, zmax, 3, 3, nz, log);
  proc::set_region(st, "silicon", -1, log);
  return st;
}

// ---------------------------------------------------------------------------
// [W-8] Screen-oxide MC attenuation / STI shielding / analytic screening
// offset: FIXED (moved to tests/test_implant_materials.cpp). All three now
// PASS -- the multi-material transport wiring itself was correct; two check
// bugs were found and fixed along the way: (1) deposit() was called with
// 0.005e-4 cm (5nm) instead of the intended 50nm oxide layer in the MC/
// analytic screening checks; (2) the STI check compared concentration in a
// narrow absolute-depth band deep in two differently-shaped tails instead of
// total dose reaching silicon (see test_implant_materials.cpp for the full
// measured writeup).
// ---------------------------------------------------------------------------
// [W-7/A-7] Resist visible in saved VTU: FIXED (moved to tests/test_photo.cpp
// test_resist_visible_in_vtu). proc::save now emits a "<base>_stack.vtu"
// sidecar (via proc::save_stack) by default when a resist stack is present.
// [W-7] save_state with resist stack: FIXED (moved to tests/test_state_io.cpp
// test_errors). save_state/export_device now warn and save/export without the
// stack instead of throwing.
// ---------------------------------------------------------------------------

// ---------------------------------------------------------------------------
// [W-3] Deck parity: `etch` and `pdbset` deck commands: FIXED (moved to
// tests/test_flow.cpp test_deck_new_commands, which also covers the other
// new deck commands added alongside them: oxidize2d/sper/mechanics/refine/
// save_state/load_state/export_device/diffuse ramp=).
// ---------------------------------------------------------------------------
// [C-1] Channeling tail in the analytic implant: dual-Pearson calibration
// should keep the analytic B 40 keV pearson profile within 10x of the MC
// (channeling) profile at depth 2*Rp. Measured today at 2*Rp (~305 nm):
// analytic=4.0e9 vs MC=5.7e17 — the single-Pearson tail underestimates the
// channeling tail by ~8 orders of magnitude.
// ---------------------------------------------------------------------------
static void check_channeling_tail(std::ostream* log) {
  SimState an = make_column(0.8e-4, 80, log);
  SimState mc = make_column(0.8e-4, 80, log);
  proc::implant_gauss(an, "B", 1e14, 40.0, 0, 0, 0, false, 0, 0, 0, 0, false,
                      "pearson", log);
  proc::implant_mc(mc, "B", 1e14, 40.0, 200000, 0, 0, 3, 0, /*channeling=*/true,
                   false, 0, 0, 0, 0, false, false, log);
  const double z0 = an.mesh.bbox().hi.z;
  const double rp = si_peak_depth(an, "B", z0);  // analytic Rp (~152 nm)

  auto conc_at = [&](const SimState& s, double depth) {
    const auto& f = s.fields.at("B");
    double best = 1e300, c = 0;
    for (std::size_t i = 0; i < f.size(); ++i) {
      const double d = std::fabs((z0 - s.mesh.cell_cent[i].z) - depth);
      if (d < best) { best = d; c = f[i]; }
    }
    return c;
  };
  const double ca = conc_at(an, 2 * rp);
  const double cm = conc_at(mc, 2 * rp);
  // "within 10x": analytic must not underestimate the MC tail by >10x.
  record("C-1", "analytic channeling tail @2Rp",
         cm > 0 && ca >= 0.1 * cm,
         fmtv("analytic=%.3e vs MC=%.3e cm^-3; want analytic >= MC/10", ca, cm));
}

// ---------------------------------------------------------------------------
// [C-3] Massoud thin-film enhancement: dry oxidation in the ~10 nm regime at
// 900 C must grow >10% MORE than the pure Deal-Grove prediction. Measured
// today: grown/DG ratio = 1.0000 (exact Deal-Grove, no thin-oxide term).
// ---------------------------------------------------------------------------
static void check_massoud(std::ostream* log) {
  SimState st = make_column(0.3e-4, 12, log);
  const double t_min = 40.0;  // ~10 nm regime at 900 C dry
  const double x = proc::oxidize(st, t_min * 60.0, 900 + 273.15, false, log);
  const double x_dg = deal_grove_step(0.0, t_min, 900.0, false) * 1e-4;
  record("C-3", "Massoud thin-oxide enhancement",
         x > 1.10 * x_dg,
         fmtv("grown=%.2f nm vs Deal-Grove=%.2f nm; want grown > 1.10*DG",
              x * 1e7, x_dg * 1e7));
}

// ---------------------------------------------------------------------------
// Tier B quantitative-benchmark checks (companion to tests/test_benchmarks.cpp,
// which holds the passing Tier A asserts and Tier C informational rows).
// These are literature-anchored targets the engine currently misses by a
// large, unambiguous factor — same WILL_FAIL contract as the checks above.
// ---------------------------------------------------------------------------

// [C-2] TED numerical stability at 750-850 C: FIXED (moved to
// tests/test_ted.cpp test 9). Root cause: the 1a implicit CI/CV solve
// undershoots at implant seed spikes (deferred non-orthogonal correction
// is explicit); a negative CI fed the BIC kinetics a negative forward
// rate, creating cluster mass from nothing. Fixed by flooring CI/CV
// after the 1a solves and clamping ratio/forward/cl_old >= 0.


// [C-2] TED enhancement magnitude: time-averaged Dt enhancement
// (sigma_ted^2-sigma0^2)/(sigma_eq^2-sigma0^2) of a damage-seeded B marker
// (1e14 cm^-2, Rp 50 nm) after 60 s at 900 C must land inside a generous
// [5x, 200x] band around the classic 10-100x marker-experiment range
// (Packan & Plummer; Stolk et al. 1997, 750-810 C short anneals; enhancement
// only decreases toward 900 C). Measured today: ~706x — the engine
// overestimates the literature enhancement by roughly an order of magnitude.
static void check_ted_enhancement_band(std::ostream* log) {
  DiffuseOpts d;
  d.temp = 1173.15;  // 900 C
  d.time = 60;
  d.verbosity = 0;
  auto sigma_of = [](const SimState& s) {
    const auto& f = s.fields.at("B");
    const double ztop = s.mesh.bbox().hi.z;
    double m = 0, md = 0, md2 = 0;
    for (std::size_t i = 0; i < f.size(); ++i) {
      const double w = f[i] * s.mesh.cell_vol[i];
      const double dd = ztop - s.mesh.cell_cent[i].z;
      m += w; md += w * dd; md2 += w * dd * dd;
    }
    const double mean = md / m;
    return std::sqrt(std::max(0.0, md2 / m - mean * mean));
  };
  SimState eq = make_column(1.0e-4, 100, log);
  proc::implant_gauss(eq, "B", 1e14, 0, 0.05e-4, 0.02e-4, 0, false, 0, 0, 0, 0,
                      false, "gauss", log);
  const double s0 = sigma_of(eq);
  proc::diffuse(eq, d, log);
  const double seq = sigma_of(eq);
  SimState td = make_column(1.0e-4, 100, log);
  proc::implant_gauss(td, "B", 1e14, 0, 0.05e-4, 0.02e-4, 0, false, 0, 0, 0, 0,
                      true, "gauss", log);
  proc::diffuse_ted(td, d, log);
  const double sted = sigma_of(td);
  const double enh = (sted * sted - s0 * s0) / (seq * seq - s0 * s0);
  record("C-2", "TED enhancement in classic band",
         std::isfinite(enh) && enh >= 5.0 && enh <= 200.0,
         fmtv("Dt enhancement=%.1fx (eq sigma->%.4g cm); want 5-200x "
              "(lit 10-100x @750-810C)", enh, seq));
}

// ---------------------------------------------------------------------------
int main() {
  std::ostringstream log;
  std::printf("SProcess-parity executable specification (IMPLEMENTATION_PLAN_v2)\n");
  std::printf("-----------------------------------------------------------------\n");

  check_channeling_tail(&log);
  check_massoud(&log);
  // Tier B quantitative benchmarks (see tests/test_benchmarks.cpp header).
  check_ted_enhancement_band(&log);

  int npass = 0;
  std::printf("\n===================== parity summary =====================\n");
  std::printf("%-8s %-40s %s\n", "task", "check", "status");
  for (const auto& r : g_results) {
    std::printf("%-8s %-40s %s\n", r.task.c_str(), r.name.c_str(),
                r.pass ? "PASS" : "FAIL");
    if (r.pass) ++npass;
  }
  std::printf("===========================================================\n");
  std::printf("%d/%zu parity checks passing — this suite is EXPECTED to fail "
              "until IMPLEMENTATION_PLAN_v2 completes\n",
              npass, g_results.size());
  if (npass < (int)g_results.size()) {
    std::printf("(registered with WILL_FAIL TRUE: ctest reports green while "
                "gaps remain; a newly-PASSing check should be moved to the "
                "regular suites)\n");
    return 1;
  }
  std::printf("all parity checks pass — remove WILL_FAIL and retire this "
              "suite into the regular tests\n");
  return 0;
}
