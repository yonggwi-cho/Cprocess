// Golden cross-process-step scenario flows (Tier 1 of the IMPLEMENTATION_PLAN_v2
// test specification). These chain multiple proc:: steps into realistic
// front-end flows and pin down invariants that hold with the CURRENT code:
// completion without throw, dose conservation, finite fields, mesh quality,
// and mask selectivity. They intentionally do NOT assert not-yet-implemented
// SProcess-parity physics (e.g. screen-oxide range attenuation) — that lives
// in tests/test_sprocess_parity.cpp (WILL_FAIL executable specification).
//
// All numeric thresholds below were calibrated by measuring the current code
// first (measured values quoted in comments) and then adding safe margin.

#include <cmath>
#include <cstdio>
#include <sstream>
#include <string>
#include <vector>

#include "cprocess/diffusion.hpp"
#include "cprocess/process.hpp"
#include "cprocess/remesh.hpp"
#include "test_util.hpp"

using namespace cp;

// ---------------------------------------------------------------------------
// Helpers
// ---------------------------------------------------------------------------

static std::string mat_of(const SimState& st, int ci) {
  auto it = st.region_material.find(st.mesh.cell_region[ci]);
  return it == st.region_material.end() ? std::string("silicon") : it->second;
}

static double total_mass(const SimState& st, const std::string& sym) {
  double m = 0;
  const auto& c = st.fields.at(sym);
  for (std::size_t i = 0; i < c.size(); ++i) m += c[i] * st.mesh.cell_vol[i];
  return m;
}

static void check_all_fields_finite(const SimState& st) {
  for (const auto& [sym, f] : st.fields)
    for (double v : f) {
      CHECK(std::isfinite(v));
      CHECK(v >= -1e-6);  // allow tiny negative numerical noise
    }
}

// Average concentration of `sp` inside vs outside the x-column band [x1,x2].
static double in_out_ratio(const SimState& st, const std::string& sp,
                           double x1, double x2) {
  const auto& f = st.fields.at(sp);
  double din = 0, vin = 0, dout = 0, vout = 0;
  for (std::size_t i = 0; i < f.size(); ++i) {
    const double x = st.mesh.cell_cent[i].x;
    const double q = f[i] * st.mesh.cell_vol[i];
    if (x >= x1 && x <= x2) { din += q; vin += st.mesh.cell_vol[i]; }
    else                    { dout += q; vout += st.mesh.cell_vol[i]; }
  }
  return (din / vin) / (dout / vout + 1e-300);
}

// ---------------------------------------------------------------------------
// Scenario 1: STI-like flow.
// mesh -> polygon trench etch -> conformal oxide fill -> blanket well implant
// (MC, window spanning the whole surface) -> short RTA (diffuse_ted).
// ---------------------------------------------------------------------------
static void scenario_sti() {
  std::printf("scenario 1: STI-like flow\n");
  std::ostringstream log;
  SimState st;
  proc::mesh_box(st, 0, 0.6e-4, 0, 0.3e-4, 0, 0.5e-4, 6, 3, 10, &log);
  proc::set_region(st, "silicon", -1, &log);
  const double z_surf = st.mesh.bbox().hi.z;

  // Trench: 0.2 um wide x 0.15 um deep, across the full y extent.
  std::vector<std::pair<double, double>> trench = {
      {0.2e-4, 0}, {0.4e-4, 0}, {0.4e-4, 0.3e-4}, {0.2e-4, 0.3e-4}};
  proc::etch(st, 0.15e-4, trench, "", &log);

  // Oxide fill: conformal deposit reaches down into the trench (a plain
  // vertical deposit() only ever adds cells above the old top surface, so it
  // cannot fill an interior trench — measured: 0 oxide cells below surface
  // with deposit(), 108 with deposit_conformal()).
  proc::deposit_conformal(st, "oxide", 0.2e-4, &log);
  int ox_below = 0;
  for (std::size_t i = 0; i < st.mesh.cells.size(); ++i)
    if (mat_of(st, (int)i) == "oxide" && st.mesh.cell_cent[i].z < z_surf)
      ++ox_below;
  std::printf("  oxide cells below original surface: %d\n", ox_below);
  CHECK(ox_below >= 20);  // measured 108

  // Well implant: MC B with an explicit window spanning the whole surface,
  // with damage seeding for the RTA.
  const BBox bb = st.mesh.bbox();
  McImplantStats mc = proc::implant_mc(st, "B", 5e12, 60.0, 40000, 0, 0, 21, 0,
                                       /*channeling=*/true,
                                       /*has_window=*/true, bb.lo.x, bb.hi.x,
                                       bb.lo.y, bb.hi.y, false,
                                       /*seed_damage=*/true, &log);
  CHECK(mc.deposited > 0);

  // Short RTA.
  const double m0 = total_mass(st, "B");
  DiffuseOpts o;
  o.temp = 1000 + 273.15;
  o.time = 15;
  o.verbosity = 0;
  proc::diffuse_ted(st, o, &log);
  const double m1 = total_mass(st, "B");
  std::printf("  B mass rel change through RTA: %.3e\n",
              std::fabs(m1 - m0) / m0);
  CHECK(std::fabs(m1 - m0) <= 2e-6 * m0);  // measured 6.1e-7

  check_all_fields_finite(st);
  const QualityStats q = mesh_quality(st.mesh);
  std::printf("  min_q=%.4f cells=%zu\n", q.min_q, st.mesh.cells.size());
  CHECK(q.min_q > 0.02);  // measured 0.4817

  std::printf("  STI-like flow passed\n");
}

// ---------------------------------------------------------------------------
// Scenario 2: LOCOS + poly gate + double-resist S/D.
// mesh -> nitride mask -> oxidize_2d (short LOCOS) -> poly gate deposit ->
// photo/mask window 1 -> MC As -> strip -> photo/mask window 2 -> MC B ->
// strip -> diffuse.
// ---------------------------------------------------------------------------
static void scenario_locos_gate_sd() {
  std::printf("scenario 2: LOCOS + poly gate + double-resist S/D\n");
  std::ostringstream log;
  SimState st;
  proc::mesh_box(st, 0, 0.8e-4, 0, 0.2e-4, 0, 0.4e-4, 8, 2, 8, &log);
  proc::set_region(st, "silicon", -1, &log);

  // Nitride mask over the active area (x in [0.3, 0.5] um).
  std::vector<std::pair<double, double>> nit = {
      {0.3e-4, 0}, {0.5e-4, 0}, {0.5e-4, 0.2e-4}, {0.3e-4, 0.2e-4}};
  proc::deposit(st, "nitride", 0.03e-4, 2, nit, &log);

  // Short 2-D LOCOS oxidation (must complete without throw).
  proc::oxidize_2d(st, 30.0, 1000 + 273.15, /*wet=*/true, &log);

  // Poly gate over the center of the active area.
  std::vector<std::pair<double, double>> gate = {
      {0.35e-4, 0}, {0.45e-4, 0}, {0.45e-4, 0.2e-4}, {0.35e-4, 0.2e-4}};
  proc::deposit(st, "poly", 0.05e-4, 2, gate, &log);

  // S/D 1: window x in [0, 0.3 um], MC As.
  proc::photo(st, 0.3e-4, 4, &log);
  proc::mask(st, 0, 0.3e-4, 0, 0.2e-4, &log);
  McImplantStats mcAs = proc::implant_mc(st, "As", 1e14, 40.0, 40000, 0, 0, 31,
                                         0, false, false, 0, 0, 0, 0, false,
                                         false, &log);
  CHECK(mcAs.deposited > 0);
  proc::strip(st, &log);
  CHECK(!st.has_stack);

  // S/D 2: a DIFFERENT window x in [0.5, 0.8 um], MC B.
  proc::photo(st, 0.3e-4, 4, &log);
  proc::mask(st, 0.5e-4, 0.8e-4, 0, 0.2e-4, &log);
  McImplantStats mcB = proc::implant_mc(st, "B", 1e14, 20.0, 40000, 0, 0, 32,
                                        0, false, false, 0, 0, 0, 0, false,
                                        false, &log);
  CHECK(mcB.deposited > 0);
  proc::strip(st, &log);
  CHECK(!st.has_stack);

  // Selectivity: each dopant must land mostly inside its own window.
  const double rAs0 = in_out_ratio(st, "As", 0, 0.3e-4);
  const double rB0 = in_out_ratio(st, "B", 0.5e-4, 0.8e-4);
  std::printf("  pre-diffuse avg-conc ratios: As in/out=%.1f B in/out=%.1f\n",
              rAs0, rB0);
  CHECK(rAs0 > 3.0);  // measured 17.9
  CHECK(rB0 > 3.0);   // measured 3.6

  // Final activation anneal: dose conservation + selectivity retained.
  const double mAs0 = total_mass(st, "As"), mB0 = total_mass(st, "B");
  DiffuseOpts o;
  o.temp = 900 + 273.15;
  o.time = 15;
  o.verbosity = 0;
  proc::diffuse(st, o, &log);
  const double dAs = std::fabs(total_mass(st, "As") - mAs0) / mAs0;
  const double dB = std::fabs(total_mass(st, "B") - mB0) / mB0;
  const double rAs1 = in_out_ratio(st, "As", 0, 0.3e-4);
  const double rB1 = in_out_ratio(st, "B", 0.5e-4, 0.8e-4);
  std::printf("  post-diffuse: dAs=%.2e dB=%.2e As in/out=%.1f B in/out=%.1f\n",
              dAs, dB, rAs1, rB1);
  CHECK(dAs <= 1e-4);  // measured 5.9e-6
  CHECK(dB <= 1e-4);   // measured 1.9e-6
  CHECK(rAs1 > 3.0);   // measured 17.9
  CHECK(rB1 > 3.0);    // measured 3.6

  check_all_fields_finite(st);
  std::printf("  LOCOS + gate + S/D flow passed\n");
}

// ---------------------------------------------------------------------------
// Scenario 3: screen oxide + implant-through + spike RTA.
// mesh -> oxidize (~20 nm) -> analytic B implant through it -> spike anneal
// via DiffuseOpts::temp_profile -> conservation / peak / solubility clamp.
//
// NOTE (by design): this scenario does NOT assert Rp attenuation by the
// screen oxide — the current implant transports through the oxide as if it
// were silicon. The physically-correct attenuation is asserted (and expected
// to fail) in test_sprocess_parity.cpp check [W-8].
// ---------------------------------------------------------------------------
static void scenario_screen_oxide_spike() {
  std::printf("scenario 3: screen oxide + spike RTA\n");
  std::ostringstream log;
  SimState st;
  proc::mesh_box(st, 0, 0.3e-4, 0, 0.3e-4, 0, 0.5e-4, 3, 3, 25, &log);
  proc::set_region(st, "silicon", -1, &log);

  // Grow ~20 nm dry screen oxide (measured 19.72 nm at 900 C / 80 min).
  const double x_ox = proc::oxidize(st, 80 * 60.0, 900 + 273.15, false, &log);
  std::printf("  screen oxide: %.2f nm\n", x_ox * 1e7);
  CHECK(x_ox > 15e-7 && x_ox < 25e-7);

  // Analytic B implant through the screen oxide.
  const double atoms = proc::implant_gauss(st, "B", 1e15, 20.0, 0, 0, 0, false,
                                           0, 0, 0, 0, &log);
  CHECK(atoms > 0);

  // Spike RTA: 900 C -> 1050 C -> 900 C over 30 s via temp_profile.
  const double m0 = total_mass(st, "B");
  DiffuseOpts o;
  o.time = 30;
  o.temp = 900 + 273.15;
  o.temp_profile = {{0, 900 + 273.15}, {15, 1050 + 273.15}, {30, 900 + 273.15}};
  o.verbosity = 0;
  proc::diffuse(st, o, &log);
  const double m1 = total_mass(st, "B");
  std::printf("  B mass rel change through spike: %.3e\n",
              std::fabs(m1 - m0) / m0);
  CHECK(std::fabs(m1 - m0) <= 1e-6 * m0);  // measured 3.3e-12
  CHECK_NEAR(st.last_temp, 900 + 273.15, 1e-6);

  // Peak depth below the Si surface: finite, positive, plausible.
  const double z_si = st.mesh.bbox().hi.z - 0.56 * x_ox;  // Si surface approx
  const auto& f = st.fields.at("B");
  double cmax = -1, z_peak = 0;
  for (std::size_t i = 0; i < f.size(); ++i)
    if (mat_of(st, (int)i) == "silicon" && f[i] > cmax) {
      cmax = f[i];
      z_peak = st.mesh.cell_cent[i].z;
    }
  const double depth = z_si - z_peak;
  std::printf("  Si peak: depth=%.1f nm cmax=%.3e\n", depth * 1e7, cmax);
  CHECK(std::isfinite(depth));
  CHECK(depth > 0);            // measured 60 nm
  CHECK(depth < 200e-7);

  // Solid-solubility clamp: active concentration capped below chemical peak.
  const auto act = proc::active_field(st, "B", -1, &log);
  double amax = 0;
  for (double v : act) {
    CHECK(std::isfinite(v));
    amax = std::max(amax, v);
  }
  std::printf("  active max=%.3e (chemical max=%.3e)\n", amax, cmax);
  CHECK(amax > 0);
  CHECK(amax < cmax);           // measured 6.76e19 < 1.56e20 (clamp engaged)
  CHECK(amax < 5e20);           // sanity: below any plausible B solubility cap

  check_all_fields_finite(st);
  std::printf("  screen oxide + spike RTA flow passed\n");
}

// ---------------------------------------------------------------------------
int main() {
  scenario_sti();
  scenario_locos_gate_sd();
  scenario_screen_oxide_spike();
  std::printf("\nall golden-flow tests passed\n");
  return 0;
}
