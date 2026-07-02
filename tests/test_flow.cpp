// Composite / integration tests chaining multiple proc:: features together:
// full front-end flow, deck-vs-proc:: equivalence, thread determinism, mass
// conservation, and TED lifecycle.

#ifdef _OPENMP
#include <omp.h>
#endif

#include <cmath>
#include <cstdio>
#include <fstream>
#include <sstream>
#include <vector>

#include "cprocess/deck.hpp"
#include "cprocess/diffusion.hpp"
#include "cprocess/process.hpp"
#include "test_util.hpp"

using namespace cp;

static bool file_nonempty(const std::string& path) {
  std::ifstream f(path, std::ios::binary | std::ios::ate);
  return f.is_open() && f.tellg() > 0;
}

static double total_mass(const SimState& st, const std::string& sym) {
  const auto& f = st.fields.at(sym);
  double m = 0;
  for (std::size_t i = 0; i < f.size(); ++i) m += f[i] * st.mesh.cell_vol[i];
  return m;
}

// ---------------------------------------------------------------------------
// Test A: full front-end flow, everything chained.
// ---------------------------------------------------------------------------
static void test_full_front_end_flow() {
  std::printf("test_full_front_end_flow\n");
  std::ostringstream log;
  SimState st;
  proc::mesh_box(st, 0, 1.0e-4, 0, 1.0e-4, 0, 0.6e-4, 8, 8, 6, &log);
  proc::set_region(st, "silicon", -1, &log);
  CHECK(!log.str().empty());

  proc::init(st, "B", 1e15, -1, &log);

  proc::implant_mc(st, "B", 2e12, 40.0, 20000, 0, 0, 123, 0,
                    /*channeling=*/true, false, 0, 0, 0, 0,
                    /*damage=*/true, &log);
  CHECK(st.fields.count("I") > 0);
  double i_after_first = 0;
  for (double v : st.fields.at("I")) i_after_first = std::max(i_after_first, v);
  CHECK(i_after_first > 0);

  proc::photo(st, 0.3e-4, 4, &log);
  CHECK(st.has_stack);

  // Hexagon approximating a circle of radius 0.25um centered at (0.5,0.5)um.
  const double cx = 0.5e-4, cy = 0.5e-4, r = 0.25e-4;
  std::vector<std::pair<double, double>> hex;
  for (int k = 0; k < 6; ++k) {
    const double ang = 2.0 * M_PI * k / 6.0;
    hex.push_back({cx + r * std::cos(ang), cy + r * std::sin(ang)});
  }
  proc::mask_polygon(st, hex, &log);

  proc::implant_mc(st, "P", 1e15, 30.0, 30000, 0, 0, 7, 0, false, false, 0, 0,
                    0, 0, false, &log);

  proc::strip(st, &log);
  CHECK(!st.has_stack);

  // Measure dose inside/outside the polygon right after strip (mesh
  // unchanged from photo/mask/strip, so cell centroids still line up).
  {
    // Compare average concentration (dose density) inside vs outside the
    // masked circle, not raw totals — the outside area is much larger than
    // the inside disc, so raw mass totals are not directly comparable.
    const auto& P = st.fields.at("P");
    double dose_in = 0, vol_in = 0, dose_out = 0, vol_out = 0;
    for (std::size_t i = 0; i < P.size(); ++i) {
      const double dx = st.mesh.cell_cent[i].x - cx;
      const double dy = st.mesh.cell_cent[i].y - cy;
      const double q = P[i] * st.mesh.cell_vol[i];
      if (dx * dx + dy * dy <= r * r) { dose_in += q; vol_in += st.mesh.cell_vol[i]; }
      else { dose_out += q; vol_out += st.mesh.cell_vol[i]; }
    }
    const double avg_in = dose_in / vol_in, avg_out = dose_out / vol_out;
    std::printf("  P avg conc in=%.3e out=%.3e\n", avg_in, avg_out);
    CHECK(avg_in > 3.0 * avg_out);
  }

  proc::deposit(st, "oxide", 0.05e-4, 2, {}, &log);
  CHECK(!log.str().empty());

  // Etch left-half rectangle.
  const BBox bb = st.mesh.bbox();
  std::vector<std::pair<double, double>> left_rect = {
      {bb.lo.x, bb.lo.y}, {(bb.lo.x + bb.hi.x) / 2, bb.lo.y},
      {(bb.lo.x + bb.hi.x) / 2, bb.hi.y}, {bb.lo.x, bb.hi.y}};
  proc::etch(st, 0.02e-4, left_rect, &log);

  const int zmax_patch = st.mesh.find_patch("zmax");
  CHECK(zmax_patch >= 0);
  proc::add_bc(st, "P", zmax_patch, 1e19, &log);

  double i_before_diffuse = 0;
  for (double v : st.fields.at("I")) i_before_diffuse = std::max(i_before_diffuse, v);

  DiffuseOpts opts;
  opts.time = 5 * 60.0;
  opts.temp = 950 + 273.15;
  proc::diffuse_ted(st, opts, &log);
  proc::clear_bc(st, &log);
  CHECK(st.bcs.empty());

  double i_after_diffuse = 0;
  for (double v : st.fields.at("I")) i_after_diffuse = std::max(i_after_diffuse, v);
  std::printf("  I max before=%.3e after=%.3e\n", i_before_diffuse, i_after_diffuse);
  CHECK(i_after_diffuse < i_before_diffuse);

  for (const auto& [sym, f] : st.fields)
    for (double v : f) {
      CHECK(std::isfinite(v));
      CHECK(v >= -1e-6);  // allow tiny negative numerical noise
    }

  const std::string vtu = "/tmp/test_flow_full.vtu";
  proc::save(st, vtu, &log);
  CHECK(file_nonempty(vtu));

  std::printf("  full front-end flow passed\n");
}

// ---------------------------------------------------------------------------
// Test B: deck vs proc:: equivalence.
// ---------------------------------------------------------------------------
static void test_deck_vs_proc_equivalence() {
  std::printf("test_deck_vs_proc_equivalence\n");

  SimState st_deck;
  {
    std::istringstream in(
        "mesh box xmax=0.4um ymax=0.4um zmax=0.4um nx=4 ny=4 nz=4\n"
        "region all material=silicon\n"
        "init species=P conc=1e15\n"
        "implant species=B dose=1e13 rp=0.1um drp=0.03um\n"
        "diffuse time=5min temp=1000C\n");
    std::ostringstream log;
    run_deck(in, st_deck, log);
  }

  SimState st_proc;
  {
    std::ostringstream log;
    proc::mesh_box(st_proc, 0, 0.4e-4, 0, 0.4e-4, 0, 0.4e-4, 4, 4, 4, &log);
    proc::set_region(st_proc, "silicon", -1, &log);
    proc::init(st_proc, "P", 1e15, -1, &log);
    proc::implant_gauss(st_proc, "B", 1e13, 0.0, 0.1e-4, 0.03e-4, 0, false, 0,
                        0, 0, 0, false, "gauss", &log);
    DiffuseOpts opts;
    opts.time = 5 * 60.0;
    opts.temp = 1000 + 273.15;
    proc::diffuse(st_proc, opts, &log);
  }

  CHECK(st_deck.fields.count("B") > 0 && st_proc.fields.count("B") > 0);
  const auto& Bd = st_deck.fields.at("B");
  const auto& Bp = st_proc.fields.at("B");
  CHECK(Bd.size() == Bp.size());
  for (std::size_t i = 0; i < Bd.size(); ++i) {
    const double scale = std::max(std::fabs(Bd[i]), std::fabs(Bp[i])) + 1.0;
    CHECK_NEAR(Bd[i], Bp[i], 1e-12 * scale);
  }

  std::printf("  deck vs proc equivalence passed\n");
}

// ---------------------------------------------------------------------------
// Test C: thread determinism.
// ---------------------------------------------------------------------------
#ifdef _OPENMP
static SimState run_gauss_diffuse(int nthreads) {
  omp_set_num_threads(nthreads);
  SimState st;
  std::ostringstream log;
  proc::mesh_box(st, 0, 0.4e-4, 0, 0.4e-4, 0, 0.4e-4, 5, 5, 5, &log);
  proc::set_region(st, "silicon", -1, &log);
  proc::init(st, "B", 1e15, -1, &log);
  proc::implant_gauss(st, "B", 1e13, 0.0, 0.1e-4, 0.03e-4, 0, false, 0, 0, 0,
                      0, false, "gauss", &log);
  const int zmax = st.mesh.find_patch("zmax");
  proc::add_bc(st, "B", zmax, 1e19, &log);
  DiffuseOpts opts;
  opts.time = 60;
  opts.temp = 1000 + 273.15;
  proc::diffuse(st, opts, &log);
  return st;
}

static SimState run_mc_amorphous(int nthreads, unsigned long long seed) {
  omp_set_num_threads(nthreads);
  SimState st;
  std::ostringstream log;
  proc::mesh_box(st, 0, 0.4e-4, 0, 0.4e-4, 0, 0.4e-4, 5, 5, 5, &log);
  proc::set_region(st, "silicon", -1, &log);
  proc::implant_mc(st, "B", 1e13, 40.0, 20000, 0, 0, seed, nthreads,
                    /*channeling=*/false, false, 0, 0, 0, 0, false, &log);
  return st;
}
#endif

static void test_thread_determinism() {
  std::printf("test_thread_determinism\n");
#ifdef _OPENMP
  SimState s1 = run_gauss_diffuse(1);
  SimState s4 = run_gauss_diffuse(4);
  const auto& B1 = s1.fields.at("B");
  const auto& B4 = s4.fields.at("B");
  CHECK(B1.size() == B4.size());
  for (std::size_t i = 0; i < B1.size(); ++i) {
    const double scale = std::max(std::fabs(B1[i]), std::fabs(B4[i])) + 1.0;
    CHECK_NEAR(B1[i], B4[i], 1e-9 * scale);
  }

  SimState m1 = run_mc_amorphous(1, 99);
  SimState m4 = run_mc_amorphous(4, 99);
  const auto& Bm1 = m1.fields.at("B");
  const auto& Bm4 = m4.fields.at("B");
  CHECK(Bm1.size() == Bm4.size());
  for (std::size_t i = 0; i < Bm1.size(); ++i) CHECK(Bm1[i] == Bm4[i]);

  omp_set_num_threads(1);
  std::printf("  thread determinism passed\n");
#else
  std::printf("  OpenMP not enabled; skipped\n");
#endif
}

// ---------------------------------------------------------------------------
// Test D: conservation chain.
// ---------------------------------------------------------------------------
static void test_conservation_chain() {
  std::printf("test_conservation_chain\n");
  std::ostringstream log;
  SimState st;
  proc::mesh_box(st, 0, 0.4e-4, 0, 0.4e-4, 0, 0.4e-4, 5, 5, 8, &log);
  proc::set_region(st, "silicon", -1, &log);
  proc::init(st, "B", 1e15, -1, &log);
  proc::implant_gauss(st, "B", 1e13, 0.0, 0.1e-4, 0.03e-4, 0, false, 0, 0, 0,
                      0, false, "gauss", &log);

  const double mass0 = total_mass(st, "B");

  DiffuseOpts opts;
  opts.time = 60;
  opts.temp = 1000 + 273.15;
  proc::diffuse(st, opts, &log);  // no BC => zero-flux everywhere

  const double mass1 = total_mass(st, "B");
  std::printf("  mass0=%.6e mass1=%.6e\n", mass0, mass1);
  CHECK_NEAR(mass1, mass0, 1e-6 * mass0);

  // Blanket etch: recompute expected mass from cells that become gas.
  const int nc = static_cast<int>(st.mesh.cells.size());
  const auto& Bpre = st.fields.at("B");
  const BBox bb = st.mesh.bbox();
  const double z_cut = bb.hi.z - 0.1e-4;
  double removed = 0, kept_expected = 0;
  for (int i = 0; i < nc; ++i) {
    const double q = Bpre[i] * st.mesh.cell_vol[i];
    if (st.mesh.cell_cent[i].z >= z_cut) removed += q;
    else kept_expected += q;
  }
  proc::etch(st, 0.1e-4, {}, &log);
  const double mass2 = total_mass(st, "B");
  std::printf("  kept_expected=%.6e mass2=%.6e\n", kept_expected, mass2);
  CHECK_NEAR(mass2, kept_expected, 1e-9 * (kept_expected + 1.0));

  proc::diffuse(st, opts, &log);
  const double mass3 = total_mass(st, "B");
  std::printf("  mass2=%.6e mass3=%.6e\n", mass2, mass3);
  CHECK_NEAR(mass3, mass2, 1e-6 * (mass2 + 1.0));

  std::printf("  conservation chain passed\n");
}

// ---------------------------------------------------------------------------
// Test E: TED lifecycle.
// ---------------------------------------------------------------------------
static double field_spread(const SimState& st, const std::string& sym) {
  const auto& f = st.fields.at(sym);
  double mean_z = 0, mass = 0;
  for (std::size_t i = 0; i < f.size(); ++i) {
    const double q = f[i] * st.mesh.cell_vol[i];
    mean_z += q * st.mesh.cell_cent[i].z;
    mass += q;
  }
  if (mass <= 0) return 0;
  mean_z /= mass;
  double var = 0;
  for (std::size_t i = 0; i < f.size(); ++i) {
    const double q = f[i] * st.mesh.cell_vol[i];
    const double dz = st.mesh.cell_cent[i].z - mean_z;
    var += q * dz * dz;
  }
  var /= mass;
  return std::sqrt(std::max(var, 0.0));
}

static void test_ted_lifecycle() {
  std::printf("test_ted_lifecycle\n");
  std::ostringstream log;
  SimState st;
  proc::mesh_box(st, 0, 0.4e-4, 0, 0.4e-4, 0, 0.4e-4, 5, 5, 10, &log);
  proc::set_region(st, "silicon", -1, &log);
  proc::init(st, "B", 1e15, -1, &log);
  proc::implant_gauss(st, "B", 1e13, 40.0, 0, 0, 0, false, 0, 0, 0, 0,
                      /*seed_damage=*/true, "gauss", &log);

  const double s0 = field_spread(st, "B");
  double i0 = 0;
  for (double v : st.fields.at("I")) i0 = std::max(i0, v);

  DiffuseOpts opts;
  opts.time = 30;
  opts.temp = 1000 + 273.15;
  proc::diffuse_ted(st, opts, &log);
  const double s1 = field_spread(st, "B");
  double i1 = 0;
  for (double v : st.fields.at("I")) i1 = std::max(i1, v);

  proc::diffuse_ted(st, opts, &log);
  const double s2 = field_spread(st, "B");
  double i2 = 0;
  for (double v : st.fields.at("I")) i2 = std::max(i2, v);

  std::printf("  spreads: s0=%.4g s1=%.4g s2=%.4g; Imax: i0=%.3e i1=%.3e i2=%.3e\n",
              s0, s1, s2, i0, i1, i2);
  const double inc1 = s1 - s0;
  const double inc2 = s2 - s1;
  CHECK(inc1 > 0);
  CHECK(inc2 < inc1);
  CHECK(i1 < i0);
  CHECK(i2 < i1);

  std::printf("  TED lifecycle passed\n");
}

// ---------------------------------------------------------------------------
int main() {
  test_full_front_end_flow();
  test_deck_vs_proc_equivalence();
  test_thread_determinism();
  test_conservation_chain();
  test_ted_lifecycle();
  std::printf("\nall flow tests passed\n");
  return 0;
}
