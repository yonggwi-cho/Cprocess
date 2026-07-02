#include <cmath>
#include <cstdio>
#include <stdexcept>
#include <vector>

#include "cprocess/materials.hpp"
#include "cprocess/process.hpp"
#include "test_util.hpp"

using namespace cp;

// Mass-weighted standard deviation of a 1-D profile along z.
static double profile_spread(const SimState& st, const std::vector<double>& c) {
  double m = 0, mz = 0, mz2 = 0;
  for (std::size_t i = 0; i < c.size(); ++i) {
    const double w = c[i] * st.mesh.cell_vol[i];
    const double z = st.mesh.cell_cent[i].z;
    m += w; mz += w * z; mz2 += w * z * z;
  }
  if (m <= 0) return 0;
  const double mean = mz / m;
  return std::sqrt(std::max(0.0, mz2 / m - mean * mean));
}

static double total_mass(const SimState& st, const std::vector<double>& c) {
  double m = 0;
  for (std::size_t i = 0; i < c.size(); ++i) m += c[i] * st.mesh.cell_vol[i];
  return m;
}

static void setup(SimState& st) {
  proc::mesh_box(st, 0, 0.3e-4, 0, 0.3e-4, 0, 1.0e-4, 4, 4, 40);
  proc::set_region(st, "silicon", -1);
  // Peak = 3e15/(sqrt(2pi)*2e-6) ~= 6.0e20 cm^-3 > C_ss(900C)=1.90e20.
  proc::implant_gauss(st, "As", 3e15, 0.0, 0.05e-4, 0.02e-4, 0.0,
                      false, 0, 0, 0, 0, false);
}

int main() {
  DiffuseOpts o;
  o.temp = 1173.15;  // 900 C
  o.time = 1800;      // 30 min
  o.verbosity = 0;

  // --- 1. clamp slows diffusion; mass is conserved either way. ---
  double spread_on, spread_off;
  {
    SimState st_on, st_off;
    setup(st_on);
    setup(st_off);
    const double m0_on = total_mass(st_on, st_on.fields.at("As"));
    const double m0_off = total_mass(st_off, st_off.fields.at("As"));

    DiffuseOpts oo = o; oo.activation = true;
    proc::diffuse(st_on, oo);
    DiffuseOpts of = o; of.activation = false;
    proc::diffuse(st_off, of);

    spread_on = profile_spread(st_on, st_on.fields.at("As"));
    spread_off = profile_spread(st_off, st_off.fields.at("As"));

    const double m1_on = total_mass(st_on, st_on.fields.at("As"));
    const double m1_off = total_mass(st_off, st_off.fields.at("As"));
    CHECK(std::fabs(m1_on - m0_on) < 1e-6 * m0_on);
    CHECK(std::fabs(m1_off - m0_off) < 1e-6 * m0_off);

    std::printf("spread: activation=on %.4g um, off %.4g um\n",
                spread_on * 1e4, spread_off * 1e4);
    CHECK(spread_off > 1.02 * spread_on);
  }

  // --- 2. active <= field, low concentration cells: exact equality. ---
  {
    SimState st;
    setup(st);
    DiffuseOpts oo = o; oo.activation = true;
    proc::diffuse(st, oo);
    auto act = proc::active_field(st, "As", 1173.15);
    const auto& c = st.fields.at("As");
    CHECK(act.size() == c.size());
    for (std::size_t i = 0; i < c.size(); ++i) {
      CHECK(act[i] <= c[i] + 1e-30);
      CHECK(act[i] <= 1.90001e20);
      if (c[i] < 1.8e20) CHECK(act[i] == c[i]);
    }
  }

  // --- 3. low concentration: active == field everywhere (bit exact). ---
  {
    SimState st;
    proc::mesh_box(st, 0, 0.3e-4, 0, 0.3e-4, 0, 1.0e-4, 4, 4, 40);
    proc::set_region(st, "silicon", -1);
    proc::init(st, "B", 1e15);
    auto act = proc::active_field(st, "B", 1273.15);
    const auto& c = st.fields.at("B");
    for (std::size_t i = 0; i < c.size(); ++i) CHECK(act[i] == c[i]);
  }

  // --- 4. temp_k default uses st.last_temp. ---
  {
    SimState st;
    setup(st);
    DiffuseOpts oo = o; oo.activation = true;
    proc::diffuse(st, oo);  // sets st.last_temp = 1173.15
    auto a1 = proc::active_field(st, "As");            // default temp_k=-1
    auto a2 = proc::active_field(st, "As", st.last_temp);
    CHECK(a1.size() == a2.size());
    for (std::size_t i = 0; i < a1.size(); ++i) CHECK(a1[i] == a2[i]);
  }

  // --- 5. missing field throws. ---
  {
    SimState st;
    proc::mesh_box(st, 0, 0.3e-4, 0, 0.3e-4, 0, 1.0e-4, 4, 4, 40);
    proc::set_region(st, "silicon", -1);
    bool threw = false;
    try {
      proc::active_field(st, "P");
    } catch (const std::runtime_error&) {
      threw = true;
    }
    CHECK(threw);
  }

  std::printf("activation tests passed\n");
  return 0;
}
