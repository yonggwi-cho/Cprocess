// Tests for proc::oxidize (P1-6): blanket 1D thermal oxidation flow.

#include <cmath>
#include <cstdio>
#include <sstream>

#include "cprocess/deck.hpp"
#include "cprocess/oxidation.hpp"
#include "cprocess/process.hpp"
#include "test_util.hpp"

using namespace cp;

// Common mesh: 0.2x0.2x0.5 um, z cell height 0.01 um.
static SimState make_state(std::ostream& log) {
  SimState st;
  proc::mesh_box(st, 0, 0.2e-4, 0, 0.2e-4, 0, 0.5e-4, 4, 4, 50, &log);
  proc::set_region(st, "silicon", -1, &log);
  return st;
}

static double total_mass(const SimState& st, const std::string& sym) {
  const auto& f = st.fields.at(sym);
  double m = 0;
  for (std::size_t i = 0; i < f.size(); ++i) m += f[i] * st.mesh.cell_vol[i];
  return m;
}

// ---------------------------------------------------------------------------
// Test 1: analytic agreement with deal_grove_step (dry 1000C, 60 min).
// ---------------------------------------------------------------------------
static void test_analytic_agreement() {
  std::printf("test_analytic_agreement\n");
  std::ostringstream log;
  SimState st = make_state(log);

  const double tox_cm = proc::oxidize(st, 3600.0, 1273.15, false, &log);
  const double tox_um = tox_cm * 1e4;

  const double analytic_um = deal_grove_step(0, 60.0, 1000.0, false);
  std::printf("  tox=%.6g um, analytic=%.6g um\n", tox_um, analytic_um);
  CHECK_NEAR(analytic_um, 0.0540, 0.0540 * 0.05);  // sanity check on constant
  CHECK(std::fabs(tox_um - analytic_um) <= 1e-9 * std::fabs(analytic_um));

  // Tag-measured oxide thickness: bbox top - z_si_top.
  const BBox bb = st.mesh.bbox();
  double z_si_top = -1e300;
  const int nc = static_cast<int>(st.mesh.cells.size());
  for (int ci = 0; ci < nc; ++ci) {
    auto it = st.region_material.find(st.mesh.cell_region[ci]);
    const std::string mat = it == st.region_material.end() ? "silicon" : it->second;
    if (mat != "silicon" && mat != "si") continue;
    for (int k = 0; k < 4; ++k)
      z_si_top = std::max(z_si_top, st.mesh.nodes[st.mesh.cells[ci][k]].z);
  }
  const double tag_tox_um = (bb.hi.z - z_si_top) * 1e4;
  std::printf("  tag-measured tox=%.6g um\n", tag_tox_um);
  CHECK(std::fabs(tag_tox_um - analytic_um) <= 0.20 * analytic_um);

  CHECK(log.str().find("[oxidize]") != std::string::npos);
  std::printf("  analytic agreement passed\n");
}

// ---------------------------------------------------------------------------
// Test 2: self-consistency of repeated application.
// ---------------------------------------------------------------------------
static void test_repeat_self_consistency() {
  std::printf("test_repeat_self_consistency\n");
  std::ostringstream log;

  SimState st_two_step = make_state(log);
  proc::oxidize(st_two_step, 1800.0, 1273.15, false, &log);
  const double tox_two = proc::oxidize(st_two_step, 1800.0, 1273.15, false, &log);

  SimState st_one_step = make_state(log);
  const double tox_one = proc::oxidize(st_one_step, 3600.0, 1273.15, false, &log);

  std::printf("  tox_two_step=%.6g cm, tox_one_step=%.6g cm\n", tox_two, tox_one);
  CHECK(std::fabs(tox_two - tox_one) <= 0.05 * tox_one);
  std::printf("  repeat self-consistency passed\n");
}

// ---------------------------------------------------------------------------
// Test 3: volume bookkeeping.
// ---------------------------------------------------------------------------
static void test_volume_balance() {
  std::printf("test_volume_balance\n");
  std::ostringstream log;
  SimState st = make_state(log);

  double vol_before = 0;
  for (double v : st.mesh.cell_vol) vol_before += v;

  const double x0_cm = 0.0;
  const double tox_cm = proc::oxidize(st, 3600.0, 1273.15, false, &log);
  const double dx_ox = tox_cm - x0_cm;
  const double rise = 0.56 * dx_ox;

  double vol_after = 0;
  for (double v : st.mesh.cell_vol) vol_after += v;

  const double area = 0.2e-4 * 0.2e-4;
  const double expected_dvol = rise * area;
  std::printf("  dvol=%.6e expected=%.6e\n", vol_after - vol_before, expected_dvol);
  CHECK(std::fabs((vol_after - vol_before) - expected_dvol) <= 0.05 * expected_dvol);
  std::printf("  volume balance passed\n");
}

// ---------------------------------------------------------------------------
// Test 4: dopant conservation across the Si->SiO2 frozen band.
// ---------------------------------------------------------------------------
static void test_dopant_conservation() {
  std::printf("test_dopant_conservation\n");
  std::ostringstream log;
  SimState st = make_state(log);
  proc::init(st, "B", 1e18, -1, &log);

  const double mass_before = total_mass(st, "B");
  const double tox_cm = proc::oxidize(st, 3600.0, 1273.15, false, &log);
  const double mass_after = total_mass(st, "B");

  std::printf("  mass_before=%.6e mass_after=%.6e\n", mass_before, mass_after);
  CHECK(std::fabs(mass_after - mass_before) <= 0.05 * mass_before);

  const double dx_ox = tox_cm;  // x0 was 0
  const double dx_si = 0.44 * dx_ox;
  const double z_si_top_orig = 0.5e-4;
  const double z_if = z_si_top_orig - dx_si;

  const auto& B = st.fields.at("B");
  const int nc = static_cast<int>(st.mesh.cells.size());
  bool above_zero = true;
  double max_in_band = 0;
  for (int ci = 0; ci < nc; ++ci) {
    const double z = st.mesh.cell_cent[ci].z;
    if (z > z_si_top_orig) {
      if (B[ci] != 0.0) above_zero = false;
    } else if (z > z_if && z <= z_si_top_orig) {
      max_in_band = std::max(max_in_band, B[ci]);
    }
  }
  CHECK(above_zero);
  std::printf("  max_in_band=%.6e\n", max_in_band);
  CHECK(max_in_band > 0.5e18);
  std::printf("  dopant conservation passed\n");
}

// ---------------------------------------------------------------------------
// Test 5: oxidize on a non-Si/SiO2 top surface throws.
// ---------------------------------------------------------------------------
static void test_gas_surface_throws() {
  std::printf("test_gas_surface_throws\n");
  std::ostringstream log;
  SimState st = make_state(log);
  // Blanket etch (P1-7) now physically removes cells, so the top surface
  // would simply be bare Si again; use a polygon etch (still gas-retag) to
  // leave a non-Si/SiO2 top surface for oxidize's guard to reject.
  const std::vector<std::pair<double, double>> full_poly = {
      {0, 0}, {0.2e-4, 0}, {0.2e-4, 0.2e-4}, {0, 0.2e-4}};
  proc::etch(st, 0.05e-4, full_poly, "", &log);

  bool threw = false;
  try {
    proc::oxidize(st, 600.0, 1273.15, false, &log);
  } catch (const std::runtime_error&) {
    threw = true;
  }
  CHECK(threw);
  std::printf("  gas surface throws passed\n");
}

// ---------------------------------------------------------------------------
// Test 6: deck command.
// ---------------------------------------------------------------------------
static void test_deck_oxidize() {
  std::printf("test_deck_oxidize\n");
  SimState st;
  std::ostringstream log;
  std::istringstream in(
      "mesh box xmax=0.2um ymax=0.2um zmax=0.5um nx=4 ny=4 nz=50\n"
      "region all material=silicon\n"
      "oxidize time=60min temp=1000C ambient=dry\n");
  run_deck(in, st, log);
  CHECK(log.str().find("[oxidize]") != std::string::npos);
  std::printf("  deck oxidize passed\n");
}

// ---------------------------------------------------------------------------
int main() {
  test_analytic_agreement();
  test_repeat_self_consistency();
  test_volume_balance();
  test_dopant_conservation();
  test_gas_surface_throws();
  test_deck_oxidize();
  std::printf("\nall oxidize flow tests passed\n");
  return 0;
}
