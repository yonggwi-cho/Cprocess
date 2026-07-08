// End-to-end integration test for the proc:: native API.
//
// Exercises mesh_box → set_region → init → implant_gauss → implant_mc →
// photo → mask → strip → diffuse → save, Dirichlet BCs, multi-species, and
// VTU output.

#include <cstdio>
#include <cstring>
#include <fstream>
#include <sstream>

#include "cprocess/deck.hpp"
#include "cprocess/diffusion.hpp"
#include "cprocess/process.hpp"
#include "test_util.hpp"

using namespace cp;

// ---------------------------------------------------------------------------
// Helpers
// ---------------------------------------------------------------------------

static bool file_nonempty(const std::string& path) {
  std::ifstream f(path, std::ios::binary | std::ios::ate);
  return f.is_open() && f.tellg() > 0;
}

// Build a basic silicon SimState: 0.4×0.4×0.8 µm, 4×4×8 cells.
static SimState make_basic_state(std::ostream& log) {
  SimState st;
  proc::mesh_box(st, 0, 0.4e-4, 0, 0.4e-4, 0, 0.8e-4, 4, 4, 8, &log);
  proc::set_region(st, "silicon", -1, &log);
  return st;
}

// ---------------------------------------------------------------------------
// Test 1: full flow — implant_gauss, implant_mc with and without a window,
//          photo/mask/strip deck, diffuse, save.
// ---------------------------------------------------------------------------
static void test_full_flow() {
  std::printf("test_full_flow\n");
  std::ostringstream log;
  SimState st = make_basic_state(log);

  CHECK(st.has_mesh);
  CHECK(static_cast<int>(st.mesh.cells.size()) > 0);

  // init background B.
  proc::init(st, "B", 1e15, -1, &log);
  CHECK(st.fields.count("B") > 0);

  // Analytic Gaussian P implant (full surface).
  const double atoms = proc::implant_gauss(st, "P", 1e13, 50.0, 0, 0, 0,
                                           false, 0, 0, 0, 0, &log);
  CHECK(atoms > 0);
  CHECK(st.fields.count("P") > 0);

  // Analytic Gaussian B implant with a lateral window (x in [0, 0.2 µm]).
  const double atoms_win = proc::implant_gauss(st, "B", 5e12, 30.0, 0, 0, 0,
                                               true, 0, 0.2e-4, 0, 0.4e-4, &log);
  CHECK(atoms_win > 0);

  // MC implant without window.
  McImplantStats mc1 = proc::implant_mc(st, "As", 2e13, 80.0, 20000,
                                         0, 0, 42, 0, false,
                                         false, 0, 0, 0, 0, false, false, &log);
  CHECK(mc1.deposited > 0);

  // MC implant with window (x in [0, 0.2 µm]).
  McImplantStats mc2 = proc::implant_mc(st, "P", 5e13, 30.0, 20000,
                                         0, 0, 7, 0, false,
                                         true, 0, 0.2e-4, 0, 0.4e-4, false,
                                         false, &log);
  CHECK(mc2.deposited > 0);

  // Diffuse.
  DiffuseOpts opts;
  opts.temp = 1273.15;  // 1000 °C
  opts.time = 600;      // 10 min
  proc::diffuse(st, opts, &log);

  const auto& B = st.fields.at("B");
  for (double v : B) CHECK(std::isfinite(v));

  // Save VTU.
  const std::string vtu = "/tmp/test_integration_full.vtu";
  proc::save(st, vtu, &log);
  CHECK(file_nonempty(vtu));

  std::printf("  full flow passed\n");
}

// ---------------------------------------------------------------------------
// Test 2: photo → mask → strip.
// ---------------------------------------------------------------------------
static void test_photo_flow() {
  std::printf("test_photo_flow\n");
  std::ostringstream log;
  // Use a slightly larger mesh so photo/strip is stable.
  SimState st;
  proc::mesh_box(st, 0, 1e-4, 0, 1e-4, 0, 0.5e-4, 8, 8, 4, &log);
  proc::set_region(st, "silicon", -1, &log);
  proc::init(st, "B", 1e15, -1, &log);

  proc::photo(st, 0.4e-4, 4, &log);
  CHECK(st.has_stack);

  // Open gate window at [0.35, 0.65 µm].
  proc::mask(st, 0, 0.35e-4, 0, 1e-4, &log);
  proc::mask(st, 0.65e-4, 1e-4, 0, 1e-4, &log);

  McImplantStats mc = proc::implant_mc(st, "P", 5e14, 30.0, 40000,
                                        0, 0, 17, 0, false,
                                        false, 0, 0, 0, 0, false, false, &log);
  CHECK(mc.deposited > 0);

  proc::strip(st, &log);
  CHECK(!st.has_stack);
  CHECK(st.fields.count("P") > 0);

  // Dose check: open region must dominate.
  const auto& P = st.fields.at("P");
  const int nc = static_cast<int>(st.mesh.cells.size());
  double dose_open = 0, dose_gate = 0;
  for (int i = 0; i < nc; ++i) {
    const double x = st.mesh.cell_cent[i].x;
    const double q = P[i] * st.mesh.cell_vol[i];
    if (x < 0.35e-4 || x > 0.65e-4) dose_open += q;
    else                             dose_gate += q;
  }
  std::printf("  P: open=%.3e gate=%.3e ratio=%.1fx\n",
              dose_open, dose_gate, dose_open / (dose_gate + 1e-30));
  CHECK(dose_open > 3.0 * dose_gate);

  std::printf("  photo flow passed\n");
}

// ---------------------------------------------------------------------------
// Test 3: Dirichlet BCs — add, diffuse, clear, diffuse again without error.
// ---------------------------------------------------------------------------
static void test_dirichlet_bc() {
  std::printf("test_dirichlet_bc\n");
  std::ostringstream log;
  SimState st = make_basic_state(log);
  proc::init(st, "P", 1e15, -1, &log);

  // Add a BC on zmax (patch 5).
  proc::add_bc(st, "P", 5, 1e20, &log);
  CHECK(st.bcs.size() == 1);

  DiffuseOpts opts;
  opts.temp = 1273.15;
  opts.time = 300;
  proc::diffuse(st, opts, &log);

  const auto& P = st.fields.at("P");
  for (double v : P) CHECK(std::isfinite(v));

  // Clear BCs and diffuse again — must not throw.
  proc::clear_bc(st, &log);
  CHECK(st.bcs.empty());

  opts.time = 60;
  bool ok = true;
  try { proc::diffuse(st, opts, &log); }
  catch (...) { ok = false; }
  CHECK(ok);

  std::printf("  dirichlet BC passed\n");
}

// ---------------------------------------------------------------------------
// Test 4: multi-species — B and As both diffuse; both fields stay finite.
// ---------------------------------------------------------------------------
static void test_multi_species() {
  std::printf("test_multi_species\n");
  std::ostringstream log;
  SimState st = make_basic_state(log);

  proc::init(st, "B",  1e15, -1, &log);
  proc::init(st, "As", 5e14, -1, &log);

  // Add a small MC implant so each species has a non-trivial profile.
  proc::implant_gauss(st, "B",  1e13, 30.0, 0, 0, 0, false, 0, 0, 0, 0, &log);
  proc::implant_gauss(st, "As", 2e13, 60.0, 0, 0, 0, false, 0, 0, 0, 0, &log);

  DiffuseOpts opts;
  opts.temp = 1273.15;
  opts.time = 300;
  proc::diffuse(st, opts, &log);

  CHECK(st.fields.count("B")  > 0);
  CHECK(st.fields.count("As") > 0);
  for (double v : st.fields.at("B"))  CHECK(std::isfinite(v));
  for (double v : st.fields.at("As")) CHECK(std::isfinite(v));

  std::printf("  multi-species passed\n");
}

// ---------------------------------------------------------------------------
// Test 5: proc_save writes a non-empty VTU file.
// ---------------------------------------------------------------------------
static void test_save_vtu() {
  std::printf("test_save_vtu\n");
  std::ostringstream log;
  SimState st = make_basic_state(log);
  proc::init(st, "B", 1e15, -1, &log);

  const std::string path = "/tmp/test_integration_save.vtu";
  proc::save(st, path, &log);
  CHECK(file_nonempty(path));

  std::printf("  proc_save VTU passed\n");
}

// ---------------------------------------------------------------------------
int main() {
  test_full_flow();
  test_photo_flow();
  test_dirichlet_bc();
  test_multi_species();
  test_save_vtu();
  std::printf("\nall integration tests passed\n");
  return 0;
}
