#pragma once
#include <cstdint>
#include <vector>

#include "materials.hpp"
#include "mesh.hpp"

namespace cp {

// One constituent element of a compound target, e.g. Si or O within SiO2.
// `x` is the number fraction (sum of x over all components of a material
// must equal 1).
struct TargetComponent {
  int z = 14;         // atomic number
  double m = 28.086;   // atomic mass [amu]
  double x = 1.0;      // number fraction
};

// A stopping target for MC ion transport. Built-in helpers below provide
// silicon (crystalline, channeling-capable), photoresist, and SiO2.
//
// True compound BCA: `comp`, when non-empty, lists the target's constituent
// elements {Z_i, M_i, x_i}. Each nuclear collision picks its partner element
// with probability proportional to x_i (see walk_ion). Electronic stopping
// still uses the Bragg rule (weighted sum over components) and the free
// flight is drawn from the total atomic density `n`.
//
// `comp.empty()` is the single-element degenerate case, equivalent to
// `comp = {{z, m, 1.0}}`, and is guaranteed bit-identical to the legacy
// effective-single-element implementation (no extra RNG draws — see
// walk_ion). `z`/`m` remain the Bragg-rule effective values for logging and
// for materials that stay single-element (photoresist, vacuum).
struct TargetMaterial {
  const char* name = "Si";
  int z = 14;            // effective atomic number
  double m = 28.086;     // effective atomic mass [amu]
  double n = 4.99e22;    // TOTAL atomic density [cm^-3] (all components)
  bool crystal_si = false;  // enables crystal channeling in this material
  std::vector<TargetComponent> comp;  // empty => single element {z, m, 1.0}
};

// Crystalline silicon substrate (channeling on when McImplantParams.channeling).
TargetMaterial target_silicon();
// Organic photoresist (DNQ-novolac), ~1.2 g/cm^3, carbon-dominated.
TargetMaterial target_photoresist();
// Thermal SiO2, ~2.2 g/cm^3. True compound BCA: comp = {Si 1/3, O 2/3}.
TargetMaterial target_oxide();
// Si3N4, ~3.1 g/cm^3. True compound BCA: comp = {Si 3/7, N 4/7}.
TargetMaterial target_nitride();
// Near-vacuum / ambient: ~1000x less dense than a solid, so ions cross it
// essentially undeflected. Use it to fill the developed (open) regions above a
// patterned resist so ions reach the true silicon surface at the right depth.
TargetMaterial target_vacuum();

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
  // Opt-in: wrap lateral exits periodically even in window mode. Default off
  // preserves legacy out_of_domain behaviour bit-for-bit. Blanket-beam mode
  // (has_window=false) already wraps unconditionally regardless of this flag.
  bool lateral_wrap = false;
  long long ions = 100000;
  int threads = 0;          // 0 = hardware concurrency
  std::uint64_t seed = 1;
  bool channeling = true;   // crystal channeling + damage accumulation (Si)

  // Physical (multi-material) masking. When cell_material is non-null, the ion
  // is tracked through whatever material occupies each cell — e.g. a patterned
  // photoresist layer that slows and stops ions before they reach the silicon,
  // including lateral straggle under the mask edge. cell_material[ci] indexes
  // into material_table; material_table[0] is the substrate. When cell_material
  // is null, the whole domain is the single crystal-Si target as before.
  const std::vector<int>* cell_material = nullptr;
  std::vector<TargetMaterial> material_table;  // [0] = substrate
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
