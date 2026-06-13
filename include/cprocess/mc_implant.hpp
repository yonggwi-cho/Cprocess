#pragma once
#include <cstdint>
#include <vector>

#include "materials.hpp"
#include "mesh.hpp"

namespace cp {

// Monte Carlo ion implantation in the binary collision approximation
// (TRIM-style, amorphous target):
//
//  - free flight of one interatomic distance L = N^(-1/3) between
//    collisions; impact parameter p = p_max sqrt(U) with pi p_max^2 L N = 1
//  - nuclear scattering from the ZBL universal potential: the exact
//    classical scattering integral, precomputed once on a log-log
//    (eps, b) table and interpolated bilinearly (corteo-style); validated
//    in the test suite against an independent quadrature
//  - Lindhard-Scharff electronic stopping  S_e = k sqrt(E)
//  - every region uses silicon stopping powers (masks are not yet
//    distinguished); the target is amorphous, so channeling is absent
//
// Parallel design: ion histories are partitioned into fixed-size chunks,
// each with its own counter-seeded RNG stream; workers accumulate integer
// per-cell hit counts and per-chunk moment partials, which are reduced
// deterministically. Results are bit-identical for any thread count.
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
McImplantStats apply_mc_implant(const Mesh& mesh,
                                const std::vector<char>& silicon_mask,
                                const McImplantParams& p,
                                std::vector<double>& conc);

namespace mc {
// Exposed for validation tests.
double sin2_half_theta(double eps, double b);  // tabulated, sin^2(theta_cm/2)
double reduced_energy_per_ev(int z1, double m1, int z2, double m2);
double screening_length_cm(int z1, int z2);
}  // namespace mc

}  // namespace cp
