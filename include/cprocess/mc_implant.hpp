#pragma once
#include <cstdint>
#include <vector>

#include "materials.hpp"
#include "mesh.hpp"

namespace cp {

// Monte Carlo ion implantation in the binary collision approximation
// (TRIM-style):
//
//  - free flight L = N^(-1/3) between collisions; impact parameter drawn
//    uniformly in [b_min, pmax] area, pi*pmax^2*L*N = 1
//  - nuclear scattering: ZBL potential, exact classical scattering integral
//    on a 160x160 log-log (eps,b) table, bilinear interpolation (corteo-style)
//  - Lindhard-Scharff electronic stopping S_e = k*sqrt(E)
//
// Crystal channeling (enabled by default):
//  - Lindhard-Robinson continuum string potential for Si <100>/<110>/<111>
//  - Channeling criterion per step: E*sin^2(psi) < U_max*(1-f_amor)
//  - f_amor = damage_density / kNamorph from Kinchin-Pease displacement counts
//  - All threads share the damage array via OpenMP atomic updates so every
//    ion sees the accumulated crystal damage in real time; this improves
//    physical accuracy at the cost of non-determinism across thread counts
//  - Amorphous mode (channeling=false): results are bit-identical for any
//    thread count (deterministic chunk-based RNG streams)
struct McImplantParams {
  const Dopant* dopant = nullptr;
  double dose = 0;          // cm^-2
  double energy_kev = 0;
  double tilt_deg = 0;      // beam tilt from the surface normal
  double rotation_deg = 0;  // azimuth of the tilt
  bool has_window = false;  // source restricted to [x1,x2]x[y1,y2]
  double x1 = 0, x2 = 0, y1 = 0, y2 = 0;
  long long ions = 100000;
  int threads = 0;          // 0 = hardware concurrency
  std::uint64_t seed = 1;
  bool channeling = true;   // crystal channeling + damage accumulation (Si)
};

struct McImplantStats {
  long long deposited = 0;      // came to rest in a silicon cell
  long long backscattered = 0;  // left through the top surface
  long long transmitted = 0;    // left through the bottom
  long long out_of_domain = 0;  // left laterally (window mode only)
  long long in_mask = 0;        // rested in a non-silicon cell
  long long unbinned = 0;       // rest point not located in any cell
  double weight = 0;            // atoms represented by one simulated ion
  double rp = 0, drp = 0;       // depth moments of deposited ions [cm]
};

// Adds the implanted profile to `conc` (cm^-3) for silicon cells.
// If `damage_conc` is non-null and channeling is enabled, it is filled with
// the displaced-atom density profile [cm^-3] (Kinchin-Pease model).
McImplantStats apply_mc_implant(const Mesh& mesh,
                                const std::vector<char>& silicon_mask,
                                const McImplantParams& p,
                                std::vector<double>& conc,
                                std::vector<double>* damage_conc = nullptr);

namespace mc {
// Exposed for validation tests.
double sin2_half_theta(double eps, double b);  // tabulated, sin^2(theta_cm/2)
double reduced_energy_per_ev(int z1, double m1, int z2, double m2);
double screening_length_cm(int z1, int z2);
}  // namespace mc

}  // namespace cp
