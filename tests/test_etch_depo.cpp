// Tests for P1-7: layer_stack bookkeeping, true blanket-etch cell removal
// (with node compaction + identity field transfer), polygon-etch gas-retag,
// and etch() material selectivity.

#include <cmath>
#include <cstdio>
#include <sstream>
#include <vector>

#include "cprocess/deck.hpp"
#include "cprocess/process.hpp"
#include "cprocess/remesh.hpp"
#include "test_util.hpp"

using namespace cp;

static double total_mass(const SimState& st, const std::string& sym) {
  const auto& f = st.fields.at(sym);
  double m = 0;
  for (std::size_t i = 0; i < f.size(); ++i) m += f[i] * st.mesh.cell_vol[i];
  return m;
}

// ---------------------------------------------------------------------------
// Test 1: multi-layer deposit — layer_stack bookkeeping + silicon-overwrite
// bug regression.
// ---------------------------------------------------------------------------
static void test_multilayer_deposit() {
  std::printf("test_multilayer_deposit\n");
  std::ostringstream log;
  SimState st;
  proc::mesh_box(st, 0, 1.0e-4, 0, 1.0e-4, 0, 0.5e-4, 4, 4, 8, &log);
  proc::set_region(st, "silicon", -1, &log);

  const double z0_top = st.mesh.bbox().hi.z;
  proc::deposit(st, "oxide", 0.1e-4, 2, {}, &log);
  proc::deposit(st, "nitride", 0.1e-4, 2, {}, &log);

  const double z1_top = st.mesh.bbox().hi.z;
  CHECK_NEAR(z1_top - z0_top, 0.2e-4, 1e-9 * 0.2e-4);

  CHECK(st.layer_stack.size() == 2);
  CHECK(st.layer_stack[0].second == "nitride");
  CHECK(st.layer_stack[1].second == "oxide");

  int n_nitride = 0, n_oxide = 0, n_wrong = 0;
  for (std::size_t i = 0; i < st.mesh.cells.size(); ++i) {
    const double z = st.mesh.cell_cent[i].z;
    const std::string mat = st.region_material.count(st.mesh.cell_region[i])
        ? st.region_material.at(st.mesh.cell_region[i]) : "silicon";
    if (z > z1_top - 0.1e-4 && z < z1_top) {
      if (mat == "nitride") ++n_nitride; else ++n_wrong;
    } else if (z > z1_top - 0.2e-4 && z < z1_top - 0.1e-4) {
      if (mat == "oxide") ++n_oxide; else ++n_wrong;
    }
  }
  CHECK(n_nitride > 0);
  CHECK(n_oxide > 0);
  // Regression: oxide layer must still read "oxide", not clobbered to
  // "silicon" by the second deposit's blanket region_material write.
  CHECK(n_wrong == 0);

  std::printf("  ok: layer_stack=[nitride,oxide], oxide tag not clobbered\n");
}

// ---------------------------------------------------------------------------
// Test 2: blanket etch = true cell removal.
// ---------------------------------------------------------------------------
static void test_blanket_etch_removes_cells() {
  std::printf("test_blanket_etch_removes_cells\n");
  std::ostringstream log;
  SimState st;
  proc::mesh_box(st, 0, 1.0e-4, 0, 1.0e-4, 0, 0.5e-4, 4, 4, 8, &log);
  proc::set_region(st, "silicon", -1, &log);
  proc::deposit(st, "oxide", 0.1e-4, 2, {}, &log);
  proc::deposit(st, "nitride", 0.1e-4, 2, {}, &log);

  const int n_before = static_cast<int>(st.mesh.cells.size());
  const double top_before = st.mesh.bbox().hi.z;
  const double cell_h = 0.1e-4 / 2;  // nz_add=2 layers per 0.1um film

  proc::etch(st, 0.1e-4, {}, "", &log);

  const int n_after = static_cast<int>(st.mesh.cells.size());
  CHECK(n_after < n_before);

  const double top_after = st.mesh.bbox().hi.z;
  CHECK(std::fabs(top_after - (top_before - 0.1e-4)) <= cell_h + 1e-12);

  CHECK(st.layer_stack.size() == 1);
  CHECK(st.layer_stack[0].second == "oxide");

  const QualityStats qs = mesh_quality(st.mesh);
  CHECK(qs.min_q > 0);
  CHECK(st.mesh.cell_vol.size() == st.mesh.cells.size());

  std::printf("  ok: n_cells %d -> %d, top %.4g -> %.4g um\n", n_before,
              n_after, top_before * 1e4, top_after * 1e4);
}

// ---------------------------------------------------------------------------
// Test 3: dopant mass conservation through blanket etch (identity transfer).
// ---------------------------------------------------------------------------
static void test_blanket_etch_mass_conservation() {
  std::printf("test_blanket_etch_mass_conservation\n");
  std::ostringstream log;
  SimState st;
  proc::mesh_box(st, 0, 0.4e-4, 0, 0.4e-4, 0, 0.4e-4, 5, 5, 8, &log);
  proc::set_region(st, "silicon", -1, &log);
  proc::init(st, "B", 1e15, -1, &log);
  proc::implant_gauss(st, "B", 1e13, 50.0, 0, 0, 0, false, 0, 0, 0, 0, false,
                      "gauss", &log);

  const double z_cut = st.mesh.bbox().hi.z - 0.05e-4;
  const auto& B = st.fields.at("B");
  double kept_expected = 0;
  for (std::size_t i = 0; i < B.size(); ++i)
    if (st.mesh.cell_cent[i].z <= z_cut)
      kept_expected += B[i] * st.mesh.cell_vol[i];

  proc::etch(st, 0.05e-4, {}, "", &log);

  const double kept_actual = total_mass(st, "B");
  std::printf("  kept_expected=%.6e kept_actual=%.6e\n", kept_expected,
              kept_actual);
  CHECK_NEAR(kept_actual, kept_expected, 1e-9 * (kept_expected + 1.0));
}

// ---------------------------------------------------------------------------
// Test 4: selective etch (material=).
// ---------------------------------------------------------------------------
static void test_selective_etch() {
  std::printf("test_selective_etch\n");
  std::ostringstream log;
  SimState st;
  proc::mesh_box(st, 0, 1.0e-4, 0, 1.0e-4, 0, 0.5e-4, 4, 4, 8, &log);
  proc::set_region(st, "silicon", -1, &log);
  proc::deposit(st, "oxide", 0.1e-4, 2, {}, &log);

  const double top_before_etch = st.mesh.bbox().hi.z;
  int n_si_before = 0;
  for (std::size_t i = 0; i < st.mesh.cells.size(); ++i) {
    const std::string mat = st.region_material.count(st.mesh.cell_region[i])
        ? st.region_material.at(st.mesh.cell_region[i]) : "silicon";
    if (mat == "silicon") ++n_si_before;
  }

  // Etch depth (0.2um) exceeds the oxide film thickness (0.1um) -> all oxide
  // removed, Si cells in the depth band untouched (material mismatch).
  proc::etch(st, 0.2e-4, {}, "oxide", &log);

  int n_oxide_after = 0, n_si_after = 0;
  for (std::size_t i = 0; i < st.mesh.cells.size(); ++i) {
    const std::string mat = st.region_material.count(st.mesh.cell_region[i])
        ? st.region_material.at(st.mesh.cell_region[i]) : "silicon";
    if (mat == "oxide") ++n_oxide_after;
    if (mat == "silicon") ++n_si_after;
  }
  CHECK(n_oxide_after == 0);
  CHECK(n_si_after == n_si_before);

  const double si_top = top_before_etch - 0.1e-4;  // original Si surface
  const double top_after = st.mesh.bbox().hi.z;
  const double cell_h = 0.5e-4 / 8.0;
  CHECK(std::fabs(top_after - si_top) <= cell_h + 1e-12);

  std::printf("  ok: oxide removed, Si cells unchanged (%d), top=%.4g um\n",
              n_si_after, top_after * 1e4);
}

// ---------------------------------------------------------------------------
// Test 5: polygon etch keeps legacy gas-retag behavior (n_cells unchanged).
// ---------------------------------------------------------------------------
static void test_polygon_etch_unchanged() {
  std::printf("test_polygon_etch_unchanged\n");
  std::ostringstream log;
  SimState st;
  proc::mesh_box(st, 0, 1.0e-4, 0, 1.0e-4, 0, 0.5e-4, 4, 4, 8, &log);
  proc::set_region(st, "silicon", -1, &log);
  proc::init(st, "B", 1e15, -1, &log);

  const int n_before = static_cast<int>(st.mesh.cells.size());
  const std::vector<std::pair<double, double>> square = {
      {0.2e-4, 0.2e-4}, {0.8e-4, 0.2e-4}, {0.8e-4, 0.8e-4}, {0.2e-4, 0.8e-4}};
  proc::etch(st, 0.1e-4, square, "", &log);

  const int n_after = static_cast<int>(st.mesh.cells.size());
  CHECK(n_after == n_before);

  const double z_cut = st.mesh.bbox().hi.z - 0.1e-4;
  int gas_tag = -1;
  for (const auto& [t, m] : st.region_material)
    if (m == "gas") { gas_tag = t; break; }
  CHECK(gas_tag >= 0);

  bool found_gas = false;
  const auto& B = st.fields.at("B");
  for (std::size_t i = 0; i < st.mesh.cells.size(); ++i) {
    const Vec3& c = st.mesh.cell_cent[i];
    const bool in_poly = c.x > 0.2e-4 && c.x < 0.8e-4 && c.y > 0.2e-4 &&
                         c.y < 0.8e-4;
    if (c.z >= z_cut && in_poly) {
      CHECK(st.mesh.cell_region[i] == gas_tag);
      CHECK(B[i] == 0.0);
      found_gas = true;
    }
  }
  CHECK(found_gas);
  std::printf("  ok: n_cells unchanged (%d), polygon cells retagged gas\n",
              n_after);
}

// ---------------------------------------------------------------------------
// Test 6: resist guard — blanket etch throws while a photoresist stack is
// present.
// ---------------------------------------------------------------------------
static void test_etch_resist_guard() {
  std::printf("test_etch_resist_guard\n");
  std::ostringstream log;
  SimState st;
  proc::mesh_box(st, 0, 1.0e-4, 0, 1.0e-4, 0, 0.5e-4, 4, 4, 8, &log);
  proc::set_region(st, "silicon", -1, &log);
  proc::photo(st, 0.3e-4, 4, &log);
  CHECK(st.has_stack);

  bool threw = false;
  try {
    proc::etch(st, 0.05e-4, {}, "", &log);
  } catch (const std::runtime_error&) {
    threw = true;
  }
  CHECK(threw);
  std::printf("  ok: blanket etch after photo throws\n");
}

int main() {
  test_multilayer_deposit();
  test_blanket_etch_removes_cells();
  test_blanket_etch_mass_conservation();
  test_selective_etch();
  test_polygon_etch_unchanged();
  test_etch_resist_guard();
  std::printf("\nall etch/depo tests passed\n");
  return 0;
}
