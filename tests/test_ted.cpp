#include <cmath>
#include <cstdio>
#include <vector>

#include "cprocess/materials.hpp"
#include "cprocess/process.hpp"
#include "test_util.hpp"

using namespace cp;

// Standard deviation (spread) of a 1-D profile along z, mass-weighted.
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

// Mass-weighted mean depth (surface = ztop, depth increases downward).
static double mean_depth(const SimState& st, const std::vector<double>& c) {
  const BBox bb = st.mesh.bbox();
  double m = 0, md = 0;
  for (std::size_t i = 0; i < c.size(); ++i) {
    const double w = c[i] * st.mesh.cell_vol[i];
    const double depth = bb.hi.z - st.mesh.cell_cent[i].z;
    m += w; md += w * depth;
  }
  return (m > 0) ? md / m : 0.0;
}

// Build a shallow boron implant near the surface and return the state.
static void setup(SimState& st, bool damage) {
  proc::mesh_box(st, 0, 0.3e-4, 0, 0.3e-4, 0, 1.0e-4, 4, 4, 40);
  proc::set_region(st, "silicon", -1);
  // Shallow Gaussian B implant (Rp = 0.05 um). damage=true seeds interstitials.
  proc::implant_gauss(st, "B", 1e14, 0.0, 0.05e-4, 0.02e-4, 0.0,
                      false, 0, 0, 0, 0, damage);
}

int main() {
  const DiffuseOpts base = [] {
    DiffuseOpts o;
    o.temp = 1173.15;  // 900 C — low enough that equilibrium diffusion is slow
    o.time = 60;       // 1 minute
    o.verbosity = 0;
    return o;
  }();

  // --- 1. TED enhances diffusion vs. an equilibrium anneal of equal budget. ---
  double spread_eq, spread_ted;
  {
    SimState st;
    setup(st, /*damage=*/false);
    const double s0 = profile_spread(st, st.fields.at("B"));
    proc::diffuse(st, base);
    spread_eq = profile_spread(st, st.fields.at("B"));
    CHECK(spread_eq >= s0);  // equilibrium diffusion still spreads a little
    std::printf("equilibrium: spread %.4g -> %.4g um\n", s0 * 1e4, spread_eq * 1e4);
  }
  {
    SimState st;
    setup(st, /*damage=*/true);
    CHECK(st.fields.count("I") == 1);           // interstitials were seeded
    double Ipeak = 0;
    for (double v : st.fields.at("I")) Ipeak = std::max(Ipeak, v);
    CHECK(Ipeak > 0);
    const double s0 = profile_spread(st, st.fields.at("B"));
    proc::diffuse_ted(st, base);
    spread_ted = profile_spread(st, st.fields.at("B"));
    std::printf("ted:         spread %.4g -> %.4g um  (Ipeak=%.3g)\n",
                s0 * 1e4, spread_ted * 1e4, Ipeak);
  }
  // The whole point of TED: the damaged anneal diffuses markedly more.
  CHECK(spread_ted > 1.3 * spread_eq);
  std::printf("TED enhancement factor (spread): %.2fx\n", spread_ted / spread_eq);

  // --- 2. Interstitial excess decays: a second TED anneal adds far less. ---
  {
    SimState st;
    setup(st, /*damage=*/true);
    DiffuseOpts o = base;
    o.time = 30;
    proc::diffuse_ted(st, o);       // burns most of the interstitial excess
    const double s1 = profile_spread(st, st.fields.at("B"));
    double Iafter = 0;
    for (double v : st.fields.at("I")) Iafter = std::max(Iafter, v);
    proc::diffuse_ted(st, o);       // interstitials now depleted -> little TED
    const double s2 = profile_spread(st, st.fields.at("B"));
    const double d1 = s1;           // spread after first half
    const double d2 = s2 - s1;      // extra spread in second half
    std::printf("transient: 1st half +%.4g um, 2nd half +%.4g um, Ires=%.3g\n",
                d1 * 1e4, d2 * 1e4, Iafter);
    // Enhancement is front-loaded: the second identical anneal spreads less.
    CHECK(d2 < d1);
  }

  // --- 3. As (vacancy-dominated, fi=0.4) is less enhanced than B (fi=1.0). ---
  {
    SimState stB, stA;
    proc::mesh_box(stB, 0, 0.3e-4, 0, 0.3e-4, 0, 1.0e-4, 4, 4, 40);
    proc::set_region(stB, "silicon", -1);
    proc::implant_gauss(stB, "B", 1e14, 0.0, 0.05e-4, 0.02e-4, 0.0,
                        false, 0, 0, 0, 0, true);
    proc::mesh_box(stA, 0, 0.3e-4, 0, 0.3e-4, 0, 1.0e-4, 4, 4, 40);
    proc::set_region(stA, "silicon", -1);
    proc::implant_gauss(stA, "As", 1e14, 0.0, 0.05e-4, 0.02e-4, 0.0,
                        false, 0, 0, 0, 0, true);
    DiffuseOpts o = base;
    o.temp = 1273.15;  // 1000 C so As moves measurably
    const double b0 = profile_spread(stB, stB.fields.at("B"));
    const double a0 = profile_spread(stA, stA.fields.at("As"));
    proc::diffuse_ted(stB, o);
    proc::diffuse_ted(stA, o);
    const double db = profile_spread(stB, stB.fields.at("B")) - b0;
    const double da = profile_spread(stA, stA.fields.at("As")) - a0;
    std::printf("pair-diffusion: B enh +%.4g um, As enh +%.4g um\n",
                db * 1e4, da * 1e4);
    CHECK(db > da);  // interstitial-dominated B feels TED more than As
  }

  // --- 4. MC damage seeding: I peak is shallower than dopant Rp. ---
  {
    SimState st;
    proc::mesh_box(st, 0, 0.3e-4, 0, 0.3e-4, 0, 0.5e-4, 4, 4, 50);
    proc::set_region(st, "silicon", -1);
    proc::implant_mc(st, "B", 1e14, 50.0, 200000, 0.0, 0.0, 1, 1,
                     /*channeling=*/true, false, 0, 0, 0, 0,
                     /*seed_damage=*/true);
    CHECK(st.fields.count("I") == 1);
    const double dI = mean_depth(st, st.fields.at("I"));
    const double dB = mean_depth(st, st.fields.at("B"));
    std::printf("mc damage: mean_depth(I)=%.4g um, mean_depth(B)=%.4g um\n",
                dI * 1e4, dB * 1e4);
    CHECK(dI < dB);
  }

  // --- 5. I integral scales ~linearly with dose. ---
  {
    auto total_I = [](double dose) {
      SimState st;
      proc::mesh_box(st, 0, 0.3e-4, 0, 0.3e-4, 0, 0.5e-4, 4, 4, 50);
      proc::set_region(st, "silicon", -1);
      proc::implant_mc(st, "B", dose, 50.0, 200000, 0.0, 0.0, 1, 1,
                       /*channeling=*/true, false, 0, 0, 0, 0,
                       /*seed_damage=*/true);
      double sum = 0;
      const auto& I = st.fields.at("I");
      for (std::size_t i = 0; i < I.size(); ++i) sum += I[i] * st.mesh.cell_vol[i];
      return sum;
    };
    const double i1 = total_I(1e14);
    const double i2 = total_I(2e14);
    const double ratio = i2 / i1;
    std::printf("mc damage dose scaling: I(1e14)=%.4g, I(2e14)=%.4g, ratio=%.3f\n",
                i1, i2, ratio);
    CHECK(ratio >= 1.8 && ratio <= 2.2);
  }

  // --- 6. MC-seeded TED still enhances diffusion vs. equilibrium. ---
  {
    SimState st_eq, st_ted;
    proc::mesh_box(st_eq, 0, 0.3e-4, 0, 0.3e-4, 0, 0.5e-4, 4, 4, 50);
    proc::set_region(st_eq, "silicon", -1);
    proc::implant_mc(st_eq, "B", 1e14, 20.0, 100000, 0.0, 0.0, 1, 1,
                     /*channeling=*/true, false, 0, 0, 0, 0,
                     /*seed_damage=*/false);
    const double s0_eq = profile_spread(st_eq, st_eq.fields.at("B"));
    proc::diffuse(st_eq, base);
    const double spread_eq_mc = profile_spread(st_eq, st_eq.fields.at("B"));

    proc::mesh_box(st_ted, 0, 0.3e-4, 0, 0.3e-4, 0, 0.5e-4, 4, 4, 50);
    proc::set_region(st_ted, "silicon", -1);
    proc::implant_mc(st_ted, "B", 1e14, 20.0, 100000, 0.0, 0.0, 1, 1,
                     /*channeling=*/true, false, 0, 0, 0, 0,
                     /*seed_damage=*/true);
    const double s0_ted = profile_spread(st_ted, st_ted.fields.at("B"));
    proc::diffuse_ted(st_ted, base);
    const double spread_ted_mc = profile_spread(st_ted, st_ted.fields.at("B"));

    std::printf("mc ted: eq %.4g -> %.4g um, ted %.4g -> %.4g um\n",
                s0_eq * 1e4, spread_eq_mc * 1e4, s0_ted * 1e4, spread_ted_mc * 1e4);
    CHECK(spread_ted_mc > 1.3 * spread_eq_mc);
  }

  // --- 7. Seed cap: I never exceeds kAmorphizationDensity. ---
  {
    SimState st;
    proc::mesh_box(st, 0, 0.3e-4, 0, 0.3e-4, 0, 0.5e-4, 4, 4, 50);
    proc::set_region(st, "silicon", -1);
    proc::implant_mc(st, "B", 1e16, 50.0, 100000, 0.0, 0.0, 1, 1,
                     /*channeling=*/true, false, 0, 0, 0, 0,
                     /*seed_damage=*/true);
    double maxI = 0;
    for (double v : st.fields.at("I")) maxI = std::max(maxI, v);
    std::printf("mc damage cap: max(I)=%.4g cm^-3 (cap=%.4g)\n",
                maxI, kAmorphizationDensity);
    CHECK(maxI <= kAmorphizationDensity * (1 + 1e-12));
  }

  // --- 8. Non-channeling MC falls back to "+1" model (I still seeded). ---
  {
    SimState st;
    proc::mesh_box(st, 0, 0.3e-4, 0, 0.3e-4, 0, 0.5e-4, 4, 4, 50);
    proc::set_region(st, "silicon", -1);
    proc::implant_mc(st, "B", 1e14, 50.0, 50000, 0.0, 0.0, 1, 1,
                     /*channeling=*/false, false, 0, 0, 0, 0,
                     /*seed_damage=*/true);
    CHECK(st.fields.count("I") == 1);
    double Ipeak = 0;
    for (double v : st.fields.at("I")) Ipeak = std::max(Ipeak, v);
    std::printf("mc no-channeling fallback: Ipeak=%.4g\n", Ipeak);
    CHECK(Ipeak > 0);
  }

  std::printf("ted tests passed\n");
  return 0;
}
