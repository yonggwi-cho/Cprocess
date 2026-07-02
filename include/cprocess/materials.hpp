#pragma once
#include <array>
#include <string>
#include <vector>

namespace cp {

constexpr double kBoltzmannEv = 8.617333262e-5;  // eV/K

// Fraction of Kinchin-Pease Frenkel pairs surviving in-cascade
// recombination; the survivors seed the excess-interstitial field for TED
// (cf. "+1" model: net excess ~ dose, i.e. ~1% of total displacements).
constexpr double kFrenkelSurvival = 0.01;

// Displacement density at which crystalline Si is fully amorphized
// [cm^-3] (~12.5% of the atomic density). Also caps the seeded
// interstitial excess: cells at the cap are amorphous and their TED
// physics differs, but capping keeps the free supersaturation bounded.
constexpr double kAmorphizationDensity = 6.25e21;

enum class DopType { donor, acceptor };

// Per-cell material id (P1-9). 0/1 stay compatible with the P1-4 cell_mat
// values; 2 (the old "frozen/other" catch-all) is synonymous with kMatGas.
enum MatId : int { kMatSi = 0, kMatOxide = 1, kMatNitride = 2, kMatPoly = 3,
                    kMatGas = 4 };

// "silicon"/"si" -> kMatSi, "oxide"/"sio2" -> kMatOxide,
// "nitride"/"si3n4" -> kMatNitride, "poly"/"polysilicon" -> kMatPoly,
// "gas" and anything unrecognized -> kMatGas. Case-insensitive.
MatId material_id(const std::string& name);

// Dopant in silicon, Fair's charged point-defect (vacancy) model:
//   D = D0 + Dminus*(n/ni) + Ddminus*(n/ni)^2 + Dplus*(p/ni)   [cm^2/s]
// with each term DX = dX * exp(-eX / kT).
struct Dopant {
  std::string name;     // canonical lowercase name, e.g. "boron"
  std::string symbol;   // display symbol, e.g. "B"
  DopType type = DopType::donor;
  int z = 0;            // atomic number
  double m = 0;         // implanted isotope mass [amu]
  double d0 = 0, e0 = 0;        // neutral
  double dm = 0, em = 0;        // single negative (donors)
  double dmm = 0, emm = 0;      // double negative
  double dp = 0, ep = 0;        // single positive (acceptors)
  double ss_pre = 0, ss_e = 0;  // solid solubility Arrhenius fit [cm^-3, eV]
  // Fraction of diffusion mediated by self-interstitials (vs vacancies).
  // Governs transient enhanced diffusion: the effective diffusivity is scaled
  // by (1 - fi) + fi * (C_I / C_I*). B, P are interstitial-dominated (~1),
  // As is mixed (~0.4), Sb is vacancy-dominated (~0.1).
  double fi = 0.5;
  // Approximate projected range/moment table {energy keV, Rp cm, dRp cm,
  // gamma (skewness), beta (kurtosis)}; override with rp=/drp= in the deck
  // for accurate work.
  std::vector<std::array<double, 5>> range;

  // ── Si/SiO2 interface segregation (P1-4) ──
  // Equilibrium segregation coefficient m(T) = C_si / C_ox at the interface:
  //   m(T) = seg_m0 * exp(-seg_e / kT)
  // Interface transport coefficient (cm/s):
  //   h(T) = seg_h0 * exp(-seg_he / kT)
  double seg_m0 = 10.0, seg_e = 0.0;    // default: Si-favoring, T-independent
  double seg_h0 = 1.0e5, seg_he = 2.0;  // ~1.2e-3 cm/s at 1000 C
  // Diffusivity in SiO2 (plain Arrhenius, cm^2/s): D_ox = dox0*exp(-eox/kT)
  double dox0 = 0.0, eox = 0.0;         // 0 => immobile inside oxide

  // ── Diffusivity outside silicon (P1-9): plain Fickian D0*exp(-E/kT) ──
  double dnit0 = 0.0, enit = 0.0;    // Si3N4: 0 => perfect barrier
  double dpoly0 = 0.0, epoly = 0.0;  // poly-Si: GB-enhanced, ~10x cryst. Si
};

// Case-insensitive lookup by name or symbol ("B", "boron", ...); nullptr if
// unknown.
const Dopant* find_dopant(const std::string& name);
const std::vector<Dopant>& dopant_table();

// Intrinsic carrier density of silicon, Morin & Maita fit:
//   ni = 3.87e16 * T^1.5 * exp(-0.605 eV / kT)   [cm^-3]
double ni_si(double temp_k);

// Diffusivity at given T and normalized electron density n/ni.
double dopant_diffusivity(const Dopant& d, double temp_k, double n_over_ni);

// Approximate solid solubility (electrically active limit) [cm^-3].
double solid_solubility(const Dopant& d, double temp_k);

// Electrically active concentration: solid-solubility clamp.
//   C_act = min(C, C_ss(T));  C_ss == 0 (no fit) means "no clamp".
double active_concentration(const Dopant& d, double conc, double temp_k);

// ── Self-interstitial point-defect model (for transient enhanced diffusion) ──
// All values are order-of-magnitude literature fits for silicon; they set the
// TED time-scale and magnitude and can be tuned. Isothermal, spatially uniform.
//
//   equilibrium concentration  C_I*(T)   [cm^-3]
//   effective diffusivity      D_I(T)    [cm^2/s]
//   bulk I-V recombination     k(T)      [1/s]  (excess relaxes as exp(-k t))
double interstitial_cstar(double temp_k);
double interstitial_diffusivity(double temp_k);
double interstitial_recomb_rate(double temp_k);

// Interpolates the range table (linear in log E). Returns false if the
// dopant has no table; clamps outside the tabulated energy range.
bool implant_range(const Dopant& d, double energy_kev, double& rp, double& drp);

// Interpolates all 4 moments (Rp, dRp, gamma, beta), linear in log E, same
// clamping behavior as implant_range.
bool implant_moments(const Dopant& d, double energy_kev, double& rp,
                     double& drp, double& gamma, double& beta);

// ── Si/SiO2 interface segregation accessors (P1-4) ──
double segregation_m(const Dopant& d, double temp_k);     // seg_m0*exp(-seg_e/kT)
double segregation_h(const Dopant& d, double temp_k);     // seg_h0*exp(-seg_he/kT)
double oxide_diffusivity(const Dopant& d, double temp_k); // dox0*exp(-eox/kT); 0 if dox0<=0

// ── Material-dependent diffusivity (P1-9) ──
// D in material `mat` at temp_k, for n/ni = nni (only used for kMatSi).
double material_diffusivity(const Dopant& d, MatId mat, double temp_k,
                            double nni);

}  // namespace cp
