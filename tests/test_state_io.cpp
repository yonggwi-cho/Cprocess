// Tests for P1-11: CPRC1 binary save/load state and 1D depth profile.

#include <cmath>
#include <cstdio>
#include <sstream>
#include <string>
#include <vector>

#include "cprocess/deck.hpp"
#include "cprocess/process.hpp"
#include "test_util.hpp"

using namespace cp;

static const char* kPath = "build_test_state.cprc";

// ---------------------------------------------------------------------------
// Test 1: roundtrip bit-identical.
// ---------------------------------------------------------------------------
static void test_roundtrip_bitexact() {
  std::printf("test_roundtrip_bitexact\n");
  std::ostringstream log;
  SimState st;
  proc::mesh_box(st, 0, 1.0e-4, 0, 1.0e-4, 0, 1.0e-4, 6, 6, 12, &log);
  proc::set_region(st, "silicon", -1, &log);
  proc::init(st, "B", 1e15, -1, &log);
  proc::implant_gauss(st, "P", 1e13, 50.0, 0, 0, 0, false, 0, 0, 0, 0, false,
                      "gauss", &log);
  proc::deposit(st, "oxide", 0.05e-4, 2, {}, &log);
  int zmax_patch = st.mesh.find_patch("zmax");
  CHECK(zmax_patch >= 0);
  proc::add_bc(st, "B", zmax_patch, 1e17, &log);

  proc::save_state(st, kPath, &log);

  SimState ld;
  proc::load_state(ld, kPath, &log);

  CHECK(ld.mesh.cells.size() == st.mesh.cells.size());
  CHECK(ld.mesh.nodes.size() == st.mesh.nodes.size());
  for (std::size_t i = 0; i < st.mesh.nodes.size(); ++i) {
    CHECK(ld.mesh.nodes[i].x == st.mesh.nodes[i].x);
    CHECK(ld.mesh.nodes[i].y == st.mesh.nodes[i].y);
    CHECK(ld.mesh.nodes[i].z == st.mesh.nodes[i].z);
  }
  for (std::size_t i = 0; i < st.mesh.cells.size(); ++i) {
    CHECK(ld.mesh.cells[i] == st.mesh.cells[i]);
    CHECK(ld.mesh.cell_region[i] == st.mesh.cell_region[i]);
  }
  CHECK(ld.fields.size() == st.fields.size());
  for (const auto& [sym, conc] : st.fields) {
    const auto& ldc = ld.fields.at(sym);
    CHECK(ldc.size() == conc.size());
    for (std::size_t i = 0; i < conc.size(); ++i) CHECK(ldc[i] == conc[i]);
  }
  CHECK(ld.region_material == st.region_material);
  CHECK(ld.mesh.region_names == st.mesh.region_names);
  CHECK(ld.mesh.patch_names == st.mesh.patch_names);
  CHECK(ld.bcs.size() == st.bcs.size());
  for (std::size_t i = 0; i < st.bcs.size(); ++i) {
    CHECK(ld.bcs[i].species == st.bcs[i].species);
    CHECK(ld.bcs[i].patch == st.bcs[i].patch);
    CHECK(ld.bcs[i].conc == st.bcs[i].conc);
  }
  CHECK(ld.last_temp == st.last_temp);
  CHECK(ld.has_stack == false);
  CHECK(ld.layer_stack.empty());

  std::remove(kPath);
  std::printf("  ok: roundtrip bit-identical\n");
}

// ---------------------------------------------------------------------------
// Test 2: solve determinism after reload.
// ---------------------------------------------------------------------------
static void test_reload_then_diffuse_matches() {
  std::printf("test_reload_then_diffuse_matches\n");
  std::ostringstream log;
  SimState st;
  proc::mesh_box(st, 0, 1.0e-4, 0, 1.0e-4, 0, 1.0e-4, 6, 6, 12, &log);
  proc::set_region(st, "silicon", -1, &log);
  proc::init(st, "B", 1e15, -1, &log);
  proc::implant_gauss(st, "P", 1e13, 50.0, 0, 0, 0, false, 0, 0, 0, 0, false,
                      "gauss", &log);

  proc::save_state(st, kPath, &log);
  SimState ld;
  proc::load_state(ld, kPath, &log);
  std::remove(kPath);

  DiffuseOpts opts;
  opts.temp = 1273.15;
  opts.time = 60;
  opts.verbosity = 0;

  proc::diffuse(st, opts, nullptr);
  proc::diffuse(ld, opts, nullptr);

  const auto& B0 = st.fields.at("B");
  const auto& B1 = ld.fields.at("B");
  CHECK(B0.size() == B1.size());
  for (std::size_t i = 0; i < B0.size(); ++i) {
    const double rel = std::fabs(B0[i] - B1[i]) / (std::fabs(B0[i]) + 1.0);
    CHECK(rel < 1e-12);
  }
  std::printf("  ok: diffuse after reload matches original\n");
}

// ---------------------------------------------------------------------------
// Test 3: profile1d.
// ---------------------------------------------------------------------------
static void test_profile1d() {
  std::printf("test_profile1d\n");
  std::ostringstream log;
  SimState st;
  const double zmax = 0.8e-4;
  proc::mesh_box(st, 0, 1.0e-4, 0, 1.0e-4, 0, zmax, 6, 6, 24, &log);
  proc::set_region(st, "silicon", -1, &log);
  const double rp = 0.1e-4, drp = 0.02e-4;
  proc::implant_gauss(st, "B", 1e13, 0.0, rp, drp, 0, false, 0, 0, 0, 0, false,
                      "gauss", &log);

  const double cx = 0.5e-4, cy = 0.5e-4;
  auto prof = proc::profile1d(st, "B", cx, cy);
  CHECK(!prof.empty());
  for (std::size_t i = 1; i < prof.size(); ++i) CHECK(prof[i].first >= prof[i - 1].first);

  std::size_t imax = 0;
  for (std::size_t i = 1; i < prof.size(); ++i)
    if (prof[i].second > prof[imax].second) imax = i;

  const double z_top = zmax;
  const double expect_z = z_top - rp;
  const double cell_h = zmax / 24.0;
  CHECK(std::fabs(prof[imax].first - expect_z) <= cell_h + 1e-12);
  std::printf("  ok: profile peak z=%.4g expect=%.4g (cell=%.4g)\n",
              prof[imax].first, expect_z, cell_h);
}

// ---------------------------------------------------------------------------
// Test 4: error cases.
// ---------------------------------------------------------------------------
static void test_errors() {
  std::printf("test_errors\n");
  std::ostringstream log;

  {
    SimState ld;
    bool threw = false;
    try {
      proc::load_state(ld, "no_such_file.cprc", &log);
    } catch (const std::runtime_error&) {
      threw = true;
    }
    CHECK(threw);
  }

  {
    const char* dummy = "build_test_state_bad.cprc";
    FILE* f = std::fopen(dummy, "wb");
    CHECK(f != nullptr);
    std::fwrite("XXXXX", 1, 5, f);
    std::fclose(f);
    SimState ld;
    bool threw = false;
    try {
      proc::load_state(ld, dummy, &log);
    } catch (const std::runtime_error&) {
      threw = true;
    }
    CHECK(threw);
    std::remove(dummy);
  }

  {
    SimState st;
    proc::mesh_box(st, 0, 1.0e-4, 0, 1.0e-4, 0, 1.0e-4, 4, 4, 4, &log);
    proc::set_region(st, "silicon", -1, &log);
    proc::photo(st, 0.3e-4, 4, &log);
    // W-7 (moved from test_sprocess_parity): save_state with a live resist
    // stack must not throw — it warns and saves the base state without the
    // stack.
    std::ostringstream wlog;
    bool threw = false;
    try {
      proc::save_state(st, kPath, &wlog);
    } catch (const std::runtime_error&) {
      threw = true;
    }
    CHECK(!threw);
    CHECK(wlog.str().find("warning") != std::string::npos);
    CHECK(wlog.str().find("NOT be serialized") != std::string::npos ||
          wlog.str().find("NOT serialized") != std::string::npos);
    // The saved base state must load and match (stack members reset).
    SimState ld;
    proc::load_state(ld, kPath, &log);
    CHECK(!ld.has_stack);
    CHECK(ld.mesh.cells.size() == st.mesh.cells.size());

    proc::strip(st, &log);
    proc::save_state(st, kPath, &log);  // still fine without a stack
    std::remove(kPath);
  }

  {
    SimState st;
    proc::mesh_box(st, 0, 1.0e-4, 0, 1.0e-4, 0, 1.0e-4, 4, 4, 4, &log);
    proc::set_region(st, "silicon", -1, &log);
    proc::init(st, "B", 1e15, -1, &log);
    bool threw = false;
    try {
      proc::profile1d(st, "Xx", 0.5e-4, 0.5e-4);
    } catch (const std::runtime_error&) {
      threw = true;
    }
    CHECK(threw);

    auto out = proc::profile1d(st, "B", -10.0, -10.0);
    CHECK(out.empty());
  }

  std::printf("  ok: error cases behave as specified\n");
}

// ---------------------------------------------------------------------------
// Test 5 (P3-d): "nickel"/"nisi" region_material strings survive a roundtrip.
// ---------------------------------------------------------------------------
static void test_roundtrip_silicide_materials() {
  std::printf("test_roundtrip_silicide_materials\n");
  std::ostringstream log;
  SimState st;
  proc::mesh_box(st, 0, 0.1e-4, 0, 0.1e-4, 0, 0.5e-4, 2, 2, 250, &log);
  proc::set_region(st, "silicon", -1, &log);
  proc::deposit(st, "nickel", 0.06e-4, 2, {}, &log);
  proc::silicide(st, "nickel", 773.15, 50, &log);

  bool has_nickel = false, has_nisi = false;
  for (const auto& [tag, mat] : st.region_material) {
    if (mat == "nickel") has_nickel = true;
    if (mat == "nisi") has_nisi = true;
  }
  CHECK(has_nickel);
  CHECK(has_nisi);

  proc::save_state(st, kPath, &log);
  SimState ld;
  proc::load_state(ld, kPath, &log);

  CHECK(ld.region_material == st.region_material);
  CHECK(ld.mesh.cell_region == st.mesh.cell_region);
  bool ld_has_nickel = false, ld_has_nisi = false;
  for (const auto& [tag, mat] : ld.region_material) {
    if (mat == "nickel") ld_has_nickel = true;
    if (mat == "nisi") ld_has_nisi = true;
  }
  CHECK(ld_has_nickel);
  CHECK(ld_has_nisi);

  std::remove(kPath);
  std::printf("  ok: nickel/nisi region_material strings survive roundtrip\n");
}

int main() {
  test_roundtrip_bitexact();
  test_reload_then_diffuse_matches();
  test_profile1d();
  test_errors();
  test_roundtrip_silicide_materials();
  std::printf("\nall state_io tests passed\n");
  return 0;
}
