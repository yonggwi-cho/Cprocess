#include <cmath>
#include <cstdio>
#include <vector>

#include "cprocess/materials.hpp"
#include "cprocess/param_db.hpp"
#include "cprocess/process.hpp"
#include "test_util.hpp"

using namespace cp;

// Mass-weighted centroid-relative standard deviation of a 1-D profile along z.
static double profile_sigma(const SimState& st, const std::vector<double>& c) {
  double m = 0, mz = 0;
  for (std::size_t i = 0; i < c.size(); ++i) {
    const double w = c[i] * st.mesh.cell_vol[i];
    m += w;
    mz += w * st.mesh.cell_cent[i].z;
  }
  const double zbar = (m > 0) ? mz / m : 0.0;
  double var = 0;
  for (std::size_t i = 0; i < c.size(); ++i) {
    const double w = c[i] * st.mesh.cell_vol[i];
    const double dz = st.mesh.cell_cent[i].z - zbar;
    var += w * dz * dz;
  }
  return (m > 0) ? std::sqrt(var / m) : 0.0;
}

// Build a 6x6x24 box with a delta-like B profile in the top 2 z-layers. The
// domain height is chosen much larger than the expected diffusion length
// (so the free surface doesn't reflect the profile) but with cells much
// thinner than that length (so the finite initial box width contributes
// negligible variance next to the diffusion-broadened variance) -- otherwise
// sigma_scaled/sigma_base drifts well away from the analytic sqrt(2Dt)
// scaling regardless of concentration.
static void setup_delta(SimState& st) {
  const double H = 0.3e-4;  // 0.3 um: >>diffusion length (~0.04 um @ 1000C/600s)
  proc::mesh_box(st, 0, 0.05e-4, 0, 0.05e-4, 0, H, 6, 6, 24);
  proc::set_region(st, "silicon", -1);
  auto& f = st.fields["B"];
  f.resize(st.mesh.cells.size(), 0.0);
  const double dz = H / 24;
  for (std::size_t i = 0; i < f.size(); ++i)
    if (st.mesh.cell_cent[i].z > H - 2 * dz) f[i] = 1e15;
}

int main() {
  // --- 1. roundtrip: set/get/all/erase ---
  {
    ParamDB::instance().clear();
    ParamDB::instance().set("B.d0", 0.074);
    CHECK(ParamDB::instance().get("B.d0", -1.0) == 0.074);
    CHECK(ParamDB::instance().all().size() == 1);
    CHECK(ParamDB::instance().erase("B.d0"));
    CHECK(ParamDB::instance().get("B.d0", -1.0) == -1.0);
    ParamDB::instance().clear();
  }

  // --- 2. fallback: unset key returns fallback exactly ---
  {
    ParamDB::instance().clear();
    CHECK(ParamDB::instance().get("P.d0", 12.5) == 12.5);
    ParamDB::instance().clear();
  }

  // --- 3. diffusivity direct verification: neutral term doubling ---
  {
    ParamDB::instance().clear();
    const Dopant* d = find_dopant("B");
    CHECK(d != nullptr);
    const double T = 1273.15;
    const double D_base = dopant_diffusivity(*d, T, 1.0);
    const double kt = kBoltzmannEv * T;
    ParamDB::instance().set("B.d0", 2.0 * d->d0);
    const double D_new = dopant_diffusivity(*d, T, 1.0);
    const double expected_delta = d->d0 * std::exp(-d->e0 / kt);
    const double rel = std::fabs((D_new - D_base) - expected_delta) /
                       std::fabs(expected_delta);
    std::printf("dopant_diffusivity term doubling: D_base=%.6g D_new=%.6g "
                "delta=%.6g expected=%.6g rel=%.3g\n",
                D_base, D_new, D_new - D_base, expected_delta, rel);
    CHECK(rel < 1e-12);
    ParamDB::instance().clear();
  }

  // --- 4. profile spread scaling: B.d0 and B.dp doubled -> sigma * sqrt(2) ---
  {
    ParamDB::instance().clear();
    DiffuseOpts o;
    o.temp = 1273.15;  // 1000 C
    o.time = 600;       // 10 min
    o.field_enh = false;
    o.verbosity = 0;

    SimState st_base;
    setup_delta(st_base);
    proc::diffuse(st_base, o);
    const double sigma_base = profile_sigma(st_base, st_base.fields.at("B"));

    const Dopant* d = find_dopant("B");
    ParamDB::instance().set("B.d0", 2.0 * d->d0);
    ParamDB::instance().set("B.dp", 2.0 * d->dp);
    SimState st_scaled;
    setup_delta(st_scaled);
    proc::diffuse(st_scaled, o);
    const double sigma_scaled = profile_sigma(st_scaled, st_scaled.fields.at("B"));
    ParamDB::instance().clear();

    const double ratio = sigma_scaled / sigma_base;
    std::printf("spread scaling: sigma_base=%.6g sigma_scaled=%.6g ratio=%.4f "
                "(expect sqrt2=%.4f)\n",
                sigma_base, sigma_scaled, ratio, std::sqrt(2.0));
    CHECK(std::fabs(ratio - std::sqrt(2.0)) / std::sqrt(2.0) < 0.05);
  }

  // --- 5. unknown key is harmless: bit-identical diffuse result ---
  {
    ParamDB::instance().clear();
    DiffuseOpts o;
    o.temp = 1273.15;
    o.time = 600;
    o.field_enh = false;
    o.verbosity = 0;

    SimState st_a;
    setup_delta(st_a);
    proc::diffuse(st_a, o);

    ParamDB::instance().set("bogus.key", 1.0);
    SimState st_b;
    setup_delta(st_b);
    proc::diffuse(st_b, o);
    ParamDB::instance().clear();

    const auto& ca = st_a.fields.at("B");
    const auto& cb = st_b.fields.at("B");
    CHECK(ca.size() == cb.size());
    // Not strictly bit-identical: OpenMP reduction order in the linear solver
    // is not guaranteed identical run-to-run (~1e-15 relative, per the
    // project-wide convention of relative tol < 1e-9 for such cases). The
    // unknown key must be numerically inert well within that noise floor.
    double maxrel = 0;
    for (std::size_t i = 0; i < ca.size(); ++i)
      maxrel = std::max(maxrel,
                        std::fabs(ca[i] - cb[i]) / (std::fabs(ca[i]) + 1e-30));
    std::printf("unknown key harmless: maxrel=%.3g\n", maxrel);
    CHECK(maxrel < 1e-9);
  }

  std::printf("params tests passed\n");
  return 0;
}
