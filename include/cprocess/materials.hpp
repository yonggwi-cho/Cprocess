#pragma once
#include <array>
#include <string>
#include <vector>

namespace cp {

class ParamDB;

// ── Runtime parameter overrides (P1-10, cp::ParamDB) ──
// The functions below read every physical constant through
// ParamDB::instance().get(key, fallback) at call time, so cp::ParamDB::set()
// can override any of the following keys (raw core units: cm^2/s, eV,
// cm^-3, 1/s, dimensionless). <Sym> is Dopant::symbol ("B", "P", "As", "Sb").
// Unknown keys are silently inert.
//
//   <Sym>.d0     / <Sym>.e0     -- dopant_diffusivity neutral term
//   <Sym>.dm     / <Sym>.em     -- dopant_diffusivity single-negative term
//   <Sym>.dmm    / <Sym>.emm    -- dopant_diffusivity double-negative term
//   <Sym>.dp     / <Sym>.ep     -- dopant_diffusivity single-positive term
//   <Sym>.ss_pre / <Sym>.ss_e   -- solid_solubility
//   <Sym>.fi                    -- TED mixing fraction (run_ted)
//   <Sym>.seg_m0 / <Sym>.seg_e  -- segregation_m
//   <Sym>.seg_h0 / <Sym>.seg_he -- segregation_h
//   <Sym>.dox0   / <Sym>.eox    -- oxide_diffusivity
//   <Sym>.dnit0  / <Sym>.enit   -- material_diffusivity (kMatNitride)
//   <Sym>.dpoly0 / <Sym>.epoly  -- material_diffusivity (kMatPoly)
//   <Sym>.seg_sil_m0 / <Sym>.seg_sil_e -- segregation_m_silicide (P3-d)
//                                   (Si/silicide equilibrium ratio C_si/C_sil)
//   <Sym>.dsil0  / <Sym>.esil   -- silicide_diffusivity (P3-d, material_diffusivity
//                                   kMatSilicide)
//   silicide.nisi.b0  / silicide.nisi.eb   -- NiSi growth-law B(T) Arrhenius
//                                   (P3-d, proc::silicide), cm^2/s / eV
//   silicide.nisi.rsi / silicide.nisi.rmet -- NiSi Si/metal consumption ratios
//                                   (fraction of grown silicide thickness)
//   silicide.tisi2.b0 / silicide.tisi2.eb  -- TiSi2 growth-law B(T) Arrhenius
//   silicide.tisi2.rsi / silicide.tisi2.rmet -- TiSi2 Si/metal consumption ratios
//   I.cstar_pre  / I.cstar_e    -- interstitial_cstar
//   I.d0         / I.e0         -- interstitial_diffusivity
//   I.krec_pre   / I.krec_e     -- interstitial_recomb_rate
//   ted.smax                    -- run_ted supersaturation cap S_max
//   ted.frenkel_survival        -- multiplier on the "+1" seed_interstitials
//                                   added amount (fallback 1.0)
//   ted.k_ci                    -- C-I sink rate (P2-8), cm^3/s: excess
//                                   interstitials are removed at rate
//                                   k_ci*C_C*psi and accumulated into the
//                                   immobile "C_cl" field (run_ted only,
//                                   only when a "C" species is present).
//                                   Default 1e-19 (see src/diffusion.cpp for
//                                   the measured calibration against the
//                                   P2-8 TED-suppression acceptance test;
//                                   deviates from the task spec's 2e-21).

constexpr double kBoltzmannEv = 8.617333262e-5;  // eV/K

// P3-e: stress-diffusion/oxidation coupling constants.
constexpr double kBoltzmannErg = 1.380649e-16;  // erg/K
constexpr double kOmegaSi = 2.0e-23;            // Si atomic volume [cm^3]

// Fraction of Kinchin-Pease Frenkel pairs surviving in-cascade
// recombination; the survivors seed the excess-interstitial field for TED
// (cf. "+1" model: net excess ~ dose, i.e. ~1% of total displacements).
constexpr double kFrenkelSurvival = 0.01;

// Displacement density at which crystalline Si is fully amorphized
// [cm^-3] (~12.5% of the atomic density). Also caps the seeded
// interstitial excess: cells at the cap are amorphous and their TED
// physics differs, but capping keeps the free supersaturation bounded.
constexpr double kAmorphizationDensity = 6.25e21;

// Silicon atomic density [cm^-3] (diamond-cubic lattice, a0=5.431 A -> 8
// atoms per 1.6e-22 cm^3 cell ~= 5.0e22 cm^-3). Used by the OED interstitial
// injection source (P2-3): the growing Si/SiO2 interface consumes Si atoms
// at rate d(dx)/dt, and a fraction `oed.theta` of that atomic flux is
// assumed to inject as excess self-interstitials into the Si just below the
// interface (see proc::oxidize in src/process.cpp).
constexpr double kNSi = 5.0e22;

// neutral (P2-8): species that carries no net charge in Si (e.g. C, F, Ge).
// Excluded from the charge-neutrality nnet sum and from field-enhancement
// (see src/diffusion.cpp's nni loops in run()/run_ted()).
enum class DopType { donor, acceptor, neutral };

// Per-cell material id (P1-9). 0/1 stay compatible with the P1-4 cell_mat
// values; 2 (the old "frozen/other" catch-all) is synonymous with kMatGas.
enum MatId : int { kMatSi = 0, kMatOxide = 1, kMatNitride = 2, kMatPoly = 3,
                    kMatGas = 4,
                    kMatSilicide = 5,  // NiSi / TiSi2 (P3-d); phase kept in
                                       // region_material string ("nisi"/"tisi2")
                    kMatMetal = 6 };  // nickel / titanium (P3-d), pre-reaction
constexpr int kMatCount = 7;

// "silicon"/"si" -> kMatSi, "oxide"/"sio2" -> kMatOxide,
// "nitride"/"si3n4" -> kMatNitride, "poly"/"polysilicon" -> kMatPoly,
// "nickel"/"ni"/"titanium"/"ti" -> kMatMetal, "nisi"/"tisi2" -> kMatSilicide,
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

  // ── Silicide (P3-d): NiSi/TiSi2 common effective values ──
  // Equilibrium segregation ratio m_sil(T) = C_si / C_silicide =
  //   seg_sil_m0*exp(-seg_sil_e/kT) (< 1 -- silicide side preferred, i.e.
  //   Si-side dose loss into the silicide).
  double seg_sil_m0 = 0.3, seg_sil_e = 0.0;
  // Effective diffusivity inside the silicide [cm^2/s]:
  //   D_sil = dsil0*exp(-esil/kT). Kept nonzero by default -- 0 would trip
  //   the per-face D>0 guard in assemble() and block segregation exchange
  //   entirely (see kSegregation branch in src/diffusion.cpp).
  double dsil0 = 1.0e-3, esil = 2.0;
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

// ── Full point-defect model (P2-1): coupled I + V + {311} clusters ──
// All quantities are read from ParamDB (see keys below), evaluated at temp_k.
// This supersedes the single-psi model above inside run_ted(); the plain
// interstitial_* functions remain for backward-compatible callers.
//
//   pd.ci_star.pre / pd.ci_star.e   -- C_I* = pre*exp(-E/kT)  [cm^-3]
//   pd.cv_star.pre / pd.cv_star.e   -- C_V* = pre*exp(-E/kT)  [cm^-3]
//   pd.di.pre       / pd.di.e       -- D_I  = pre*exp(-E/kT)  [cm^2/s]
//   pd.dv.pre       / pd.dv.e       -- D_V  = pre*exp(-E/kT)  [cm^2/s]
//   pd.kbulk.factor                 -- k_bulk = factor*4*pi*a_Si*(D_I+D_V)
//   pd.c311.ktrap.factor            -- k_trap = factor*4*pi*a_Si*D_I
//   pd.c311.nu0                     -- {311} emission attempt frequency [1/s]
//   pd.c311.eb                      -- {311} binding energy [eV]
//   pd.damage.v_fraction            -- initial V/I seed split (see seed logic)
//   pd.c311.cref                    -- reference density [cm^-3] turning
//                                       k_trap [cm^3/s] into an effective
//                                       first-order {311} trap rate
//                                       k_trap_eff = k_trap*pd.c311.cref
//                                       (run_ted only; see src/diffusion.cpp)
struct PointDefectParams {
  double ci_star = 0, cv_star = 0;  // cm^-3 (equilibrium, at T)
  double d_i = 0, d_v = 0;          // cm^2/s
  double k_bulk = 0;                // cm^3/s (bimolecular I-V recombination)
  double k_trap = 0;                // cm^3/s (I capture by {311}, nominal)
  double k_emit = 0;                // 1/s ({311} emission rate = nu0*exp(-eb/kT))
};
PointDefectParams point_defect_params(double temp_k, const ParamDB& db);

// ── Dopant clustering (P2-2): immobile, electrically-inactive cluster field ──
// Single effective-composition cluster field per clustering-capable dopant:
//   B  -> BIC ("B3I", m=3 B atoms per n=1 captured interstitial, approximated
//         as a single field of clustered B atoms; forward order p=2 in
//         C_B_act, "boron pairs nucleate the cluster")
//   As -> As4V; the true 4th-order (As4) nucleation kinetics is far stiffer
//         than the reaction sub-cycle step can resolve, so it is approximated
//         as effective 2nd order in C_As_act (documented simplification, see
//         docs/tasks/P2-2_dopant_clusters.md)
// Other dopants (P, Sb) are not modeled: kf == 0 means "no clustering".
//   cl.b.kf / cl.b.nu0 / cl.b.eb   -- B forward rate [1/s] / reverse Arrhenius
//   cl.as.kf / cl.as.nu0 / cl.as.eb -- As forward rate [1/s] / reverse Arrhenius
struct ClusterParams {
  double kf = 0;        // 1/s; forward-rate prefactor (0 => not modeled)
  double kr = 0;         // 1/s; Arrhenius reverse (dissolution) rate at temp_k
  double pd_frac = 0;    // point-defect atoms exchanged per clustered dopant atom
  bool uses_v = false;   // true: forward driven by C_V/C_V* (As); false: C_I/C_I* (B)
};
ClusterParams cluster_params(const std::string& symbol, double temp_k,
                             const ParamDB& db);

// Reference concentration [cm^-3] normalizing the C_act^2 forward-rate term
// to a 1/s scale (dimensional bookkeeping constant, fixed per spec).
constexpr double kClusterCref = 1.0e20;

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

// ── Si/silicide interface segregation accessors (P3-d) ──
double segregation_m_silicide(const Dopant& d, double temp_k);  // seg_sil_m0*exp(-seg_sil_e/kT)
double silicide_diffusivity(const Dopant& d, double temp_k);    // dsil0*exp(-esil/kT)

// P3-e: per-species stress activation volume [cm^3] for the diffusion
// pressure coupling D -> D*exp(-p*V_act/kT). ParamDB key "stress.vact.<Sym>"
// overrides the built-in B/P/As defaults; any other species defaults to 0
// (no correction, exp(0)=1).
double stress_activation_volume(const Dopant& d);

// ── Material-dependent diffusivity (P1-9) ──
// D in material `mat` at temp_k, for n/ni = nni (only used for kMatSi).
double material_diffusivity(const Dopant& d, MatId mat, double temp_k,
                            double nni);

}  // namespace cp
