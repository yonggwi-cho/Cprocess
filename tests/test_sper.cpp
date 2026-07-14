// Tests for P3-c: proc::sper -- solid-phase epitaxial regrowth (SPER) of
// MC-amorphized silicon columns.

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <sstream>
#include <vector>

#include "cprocess/materials.hpp"
#include "cprocess/param_db.hpp"
#include "cprocess/process.hpp"
#include "test_util.hpp"

using namespace cp;

// 0.1 x 0.1 x 1.0 um, 2x2x200 -> 5nm cell height.
static void make_mesh(SimState& st, std::ostream& log) {
  proc::mesh_box(st, 0, 0.1e-4, 0, 0.1e-4, 0, 1.0e-4, 2, 2, 200, &log);
  proc::set_region(st, "silicon", -1, &log);
}

// Set st.fields["damage"] = dens for cells with depth (z_top - z) <= depth_cm,
// else 0.
static void set_amorphous(SimState& st, double depth_cm, double dens) {
  auto& D = st.fields["damage"];
  D.assign(st.mesh.cells.size(), 0.0);
  const double z_top = st.mesh.bbox().hi.z;
  for (std::size_t ci = 0; ci < st.mesh.cells.size(); ++ci) {
    if (z_top - st.mesh.cell_cent[ci].z <= depth_cm) D[ci] = dens;
  }
}

// Residual amorphous thickness, measured as (total amorphous volume) /
// (mesh xy footprint area). This mesh's lateral columns are all seeded
// identically (set_amorphous has no x/y dependence) and sper's per-column
// evolution is likewise laterally uniform, so this volumetric measure is
// equivalent to (and more robust than) picking a single column's z-sorted
// tet run -- individual tets within one brick layer do not share exactly
// the same centroid z, so a naive centroid sort/dedup mis-measures the
// per-layer run length.
static double residual_amorphous_thickness_cm(SimState& st, double c_am) {
  const auto& D = st.fields.at("damage");
  double vol_am = 0.0;
  for (std::size_t ci = 0; ci < st.mesh.cells.size(); ++ci)
    if (D[ci] >= c_am) vol_am += st.mesh.cell_vol[ci];
  const BBox bb = st.mesh.bbox();
  const double area = (bb.hi.x - bb.lo.x) * (bb.hi.y - bb.lo.y);
  return vol_am / area;
}

// ---------------------------------------------------------------------------
// Test 1: Arrhenius fit (v(T) at 3 temperatures).
// ---------------------------------------------------------------------------
static void test_arrhenius_fit() {
  std::printf("test_arrhenius_fit\n");
  std::ostringstream log;
  const double h = 5e-7;  // 5nm
  const double kB = kBoltzmannEv;
  const double v0 = 3.1e8, ea = 3.1;
  const double c_am = 6.25e21;

  const double temps[3] = {823.15, 873.15, 923.15};
  const double t_T[3] = {9.23e5, 7.56e4, 8.13e3};  // s, tuned for 300nm regrowth
  double v_meas[3];

  for (int k = 0; k < 3; ++k) {
    SimState st; make_mesh(st, log);
    set_amorphous(st, 0.5e-4, 1e22);  // 500nm amorphous
    const double v_analytic = v0 * std::exp(-ea / (kB * temps[k]));
    proc::sper(st, temps[k], t_T[k], &log);
    const double d_a = residual_amorphous_thickness_cm(st, c_am);
    const double regrown = 0.5e-4 - d_a;
    v_meas[k] = regrown / t_T[k];
    const double rel = std::fabs(v_meas[k] - v_analytic) / v_analytic;
    std::printf("  T=%.2f K: v_analytic=%.4g, v_meas=%.4g cm/s, rel=%.4g\n",
                temps[k], v_analytic, v_meas[k], rel);
    CHECK(rel < 0.02);
  }

  // 3-point least-squares fit of ln(v) = ln(v0) - Ea/(kB*T).
  double sx = 0, sy = 0, sxx = 0, sxy = 0;
  for (int k = 0; k < 3; ++k) {
    const double x = 1.0 / (kB * temps[k]);
    const double y = std::log(v_meas[k]);
    sx += x; sy += y; sxx += x * x; sxy += x * y;
  }
  const double n = 3.0;
  const double slope = (n * sxy - sx * sy) / (n * sxx - sx * sx);
  const double ea_fit = -slope;
  std::printf("  Ea_fit=%.6g eV (expect 3.1)\n", ea_fit);
  CHECK(std::fabs(ea_fit - ea) / ea < 0.05);
}

// ---------------------------------------------------------------------------
// Test 2: 100nm full-regrowth time window at 600 C.
// ---------------------------------------------------------------------------
static void test_full_regrowth_window() {
  std::printf("test_full_regrowth_window\n");
  std::ostringstream log;
  const double c_am = 6.25e21;
  const double temp_k = 873.15;

  // Full regrowth within 2.60e4 s.
  {
    SimState st; make_mesh(st, log);
    set_amorphous(st, 0.1e-4, 1e22);
    proc::sper(st, temp_k, 2.60e4, &log);
    double dmax = 0.0;
    for (double v : st.fields.at("damage")) dmax = std::max(dmax, v);
    std::printf("  t=2.60e4s: max damage=%.4g (expect < %.4g)\n", dmax, c_am);
    CHECK(dmax < c_am);
  }

  // Not-yet-complete at 2.40e4 s: ~4.8nm residual.
  {
    SimState st; make_mesh(st, log);
    set_amorphous(st, 0.1e-4, 1e22);
    proc::sper(st, temp_k, 2.40e4, &log);
    const double d_a = residual_amorphous_thickness_cm(st, c_am);
    const double expected = 1e-5 - 3.966e-10 * 2.40e4;  // ~4.8e-7 cm
    std::printf("  t=2.40e4s: residual=%.4g cm (expect ~%.4g cm)\n", d_a,
                expected);
    CHECK(d_a > 0);
    CHECK_NEAR(d_a, expected, 5e-7);  // within 1 cell (5nm)
  }
}

// ---------------------------------------------------------------------------
// Test 3: As metastable activation.
// ---------------------------------------------------------------------------
static void test_metastable_activation() {
  std::printf("test_metastable_activation\n");
  std::ostringstream log;
  const double temp_k = 873.15;
  const double css_as = 1.3e23 * std::exp(-0.66 / (kBoltzmannEv * temp_k));
  std::printf("  C_ss(As,600C)=%.6g\n", css_as);

  SimState st; make_mesh(st, log);
  set_amorphous(st, 0.1e-4, 1e22);
  auto& As = st.fields["As"];
  As.assign(st.mesh.cells.size(), 1e20);
  proc::sper(st, temp_k, 2.60e4, &log);  // full regrowth of the top 100nm

  const auto act = proc::active_field(st, "As", temp_k, &log);
  const auto& rg = st.fields.at("regrown");
  int n_regrown = 0, n_deep = 0;
  for (std::size_t ci = 0; ci < st.mesh.cells.size(); ++ci) {
    if (rg[ci] > 0.5) {
      ++n_regrown;
      CHECK(act[ci] > css_as);
      CHECK_NEAR(act[ci], 1e20, 1e-12 * 1e20);
    } else {
      ++n_deep;
      CHECK_NEAR(act[ci], css_as, 0.01 * css_as);
    }
  }
  std::printf("  regrown cells=%d, deep cells=%d\n", n_regrown, n_deep);
  CHECK(n_regrown > 0);
  CHECK(n_deep > 0);

  proc::set_param(st, "sper.act_factor", 1.0, &log);
  const auto act2 = proc::active_field(st, "As", temp_k, &log);
  for (std::size_t ci = 0; ci < st.mesh.cells.size(); ++ci) {
    if (rg[ci] > 0.5) CHECK_NEAR(act2[ci], css_as, 0.01 * css_as);
  }
  ParamDB::instance().clear();
  std::printf("  ok\n");
}

// ---------------------------------------------------------------------------
// Test 4: I zeroing and TED pre-processing.
// ---------------------------------------------------------------------------
static void test_i_zero_and_ted() {
  std::printf("test_i_zero_and_ted\n");
  std::ostringstream log;
  const double temp_k = 873.15;

  SimState st; make_mesh(st, log);
  set_amorphous(st, 0.1e-4, 1e22);
  auto& I = st.fields["I"];
  I.assign(st.mesh.cells.size(), 1e20);

  proc::sper(st, temp_k, 2.60e4, &log);  // full regrowth

  const auto& rg = st.fields.at("regrown");
  int n_regrown = 0, n_deep = 0;
  for (std::size_t ci = 0; ci < st.mesh.cells.size(); ++ci) {
    if (rg[ci] > 0.5) { CHECK(I[ci] == 0.0); ++n_regrown; }
    else { CHECK(I[ci] == 1e20); ++n_deep; }
  }
  std::printf("  regrown=%d (I=0), deep=%d (I unchanged)\n", n_regrown, n_deep);
  CHECK(n_regrown > 0);
  CHECK(n_deep > 0);

  auto& B = st.fields["B"];
  B.assign(st.mesh.cells.size(), 1e18);

  DiffuseOpts o;
  o.temp = 1173.15;
  o.time = 60.0;
  o.verbosity = 0;
  bool threw = false;
  try { proc::diffuse_ted(st, o, &log); }
  catch (...) { threw = true; }
  CHECK(!threw);
  std::printf("  diffuse_ted completed without exception\n");
}

// ---------------------------------------------------------------------------
// Test 5: no-op / error paths.
// ---------------------------------------------------------------------------
static void test_errors() {
  std::printf("test_errors\n");
  std::ostringstream log;

  // (a) no damage field: no-op, fields unchanged, log mentions "no amorphous".
  {
    SimState st; make_mesh(st, log);
    auto& B = st.fields["B"];
    B.assign(st.mesh.cells.size(), 1e18);
    const std::vector<double> before = B;
    std::ostringstream l2;
    proc::sper(st, 873.15, 600, &l2);
    CHECK(B == before);
    CHECK(l2.str().find("no amorphous") != std::string::npos);
  }

  // (b) time_s < 0 throws.
  {
    SimState st; make_mesh(st, log);
    set_amorphous(st, 0.1e-4, 1e22);
    bool threw = false;
    try { proc::sper(st, 873.15, -1.0, &log); }
    catch (const std::runtime_error&) { threw = true; }
    CHECK(threw);
  }

  // (c) has_stack (after photo) throws.
  {
    SimState st; make_mesh(st, log);
    set_amorphous(st, 0.1e-4, 1e22);
    proc::photo(st, 0.5e-4, 4, &log);
    bool threw = false;
    try { proc::sper(st, 873.15, 600, &log); }
    catch (const std::runtime_error&) { threw = true; }
    CHECK(threw);
    // (d) after strip: succeeds.
    proc::strip(st, &log);
    bool threw2 = false;
    try { proc::sper(st, 873.15, 600, &log); }
    catch (...) { threw2 = true; }
    CHECK(!threw2);
  }

  std::printf("  ok\n");
}

// ---------------------------------------------------------------------------
// Test 6: MC integration.
// ---------------------------------------------------------------------------
static void test_mc_integration() {
  std::printf("test_mc_integration\n");
  std::ostringstream log;
  SimState st;
  proc::mesh_box(st, 0, 0.2e-4, 0, 0.2e-4, 0, 0.5e-4, 4, 4, 50, &log);
  proc::set_region(st, "silicon", -1, &log);

  proc::implant_mc(st, "As", 1e15, 30.0, 20000, 0.0, 0.0, 12345ULL, 0, true,
                   false, 0, 0, 0, 0, false, true, &log);

  auto dit = st.fields.find("damage");
  CHECK(dit != st.fields.end());
  double dmax = 0.0;
  for (double v : dit->second) dmax = std::max(dmax, v);
  std::printf("  post-implant max damage=%.4g\n", dmax);
  CHECK(dmax > 0);

  if (dmax >= 6.25e21) {
    int n_amorphous_before = 0;
    for (double v : dit->second) if (v >= 6.25e21) ++n_amorphous_before;
    proc::sper(st, 873.15, 3.0e4, &log);
    int n_amorphous_after = 0;
    for (double v : st.fields.at("damage")) if (v >= 6.25e21) ++n_amorphous_after;
    std::printf("  amorphous cells: before=%d, after=%d\n",
                n_amorphous_before, n_amorphous_after);
    CHECK(n_amorphous_after < n_amorphous_before);
  } else {
    std::printf("  (max damage below amorphization threshold; skipping regrowth check)\n");
  }
}

int main() {
  test_arrhenius_fit();
  test_full_regrowth_window();
  test_metastable_activation();
  test_i_zero_and_ted();
  test_errors();
  test_mc_integration();
  std::printf("\nall sper tests passed\n");
  return 0;
}
