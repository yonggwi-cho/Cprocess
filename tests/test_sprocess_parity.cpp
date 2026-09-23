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
// [C-1] Channeling tail in the analytic implant: FIXED (moved to
// tests/test_dual_pearson.cpp). profile="dual" (primary Pearson-IV +
// exponential channeling-tail, calibrated against proc::implant_mc) keeps
// the analytic B 40 keV profile within 10x of the MC (channeling) profile
// at depth 2*Rp; see docs/tasks/C1_implant_moments.md.

// [C-3] Massoud thin-film enhancement: FIXED (moved to
// tests/test_oxidation.cpp). Deal-Grove B/A Arrhenius constants are now
// ParamDB-ified (ox.dry.*/ox.wet.*, bit-identical defaults) and a Massoud
// (1985) thin-oxide growth-rate enhancement term is available via
// ox.massoud.c/ox.massoud.l (default off -> bit-identical; C=0.9, L=10 nm
// gives grown/DG = 1.42 at 10 nm/900C dry, >10% target). Also added
// pressure_atm/hcl_frac/orient optional parameters to proc::oxidize()
// (all default to values reproducing current behavior exactly). See
// docs/tasks/C3_oxidation_calibration.md.
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


// [C-2] TED enhancement magnitude: FIXED (moved to tests/test_ted.cpp,
// test_ted_enhancement_classic_band). Root cause of the over-enhancement was
// investigated via several global point-defect-kinetics knobs (including an
// initial attempt exposing step_once_ted's inert 1e4 diffusivity-enhancement
// safety cap as ted.max_dv_scale, recalibrated to 500) -- that cap approach
// gave the right classic-band number in isolation but broke test_rta's
// ramp-anneal mass conservation (0.36% -> 5.5-7% error) at every cap value
// tried, so it was reverted. The adopted fix instead calibrates the B
// cluster (BIC) dissociation barrier cl.b.eb (materials.cpp), 2.7 -> 2.8 eV
// -- B-specific and confined to the clustering-consumption pathway, so it
// doesn't touch the general diffusivity-enhancement machinery RTA depends
// on. At 900C/60s the measured enhancement factor moved 298.9x -> 107.8x,
// mid-band of the check's [5,200]x window and close to the cited 10-100x
// literature range (Packan & Plummer; Stolk et al. 1997). This was the last
// remaining check in this suite (parity now 0/0 -- retired below per the
// harness's own final-summary comment: "remove WILL_FAIL and retire this
// suite into the regular tests"). See docs/tasks/README.md for the closing
// summary of the SProcess-parity section and docs/tasks/C2_ted_calibration.md
// for the full calibration writeup, including the rejected cap approach and
// the C-2 Rs/Xj extraction functions added alongside.

// ---------------------------------------------------------------------------
int main() {
  std::ostringstream log;
  std::printf("SProcess-parity executable specification (IMPLEMENTATION_PLAN_v2)\n");
  std::printf("-----------------------------------------------------------------\n");
  std::printf("(retired: every check that once lived here has been fixed and\n"
              "moved into the regular test suites -- see the per-check FIXED\n"
              "comments above and docs/tasks/README.md's closing summary.)\n");

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
