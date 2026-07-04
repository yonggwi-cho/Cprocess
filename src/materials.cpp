#include "cprocess/materials.hpp"

#include <algorithm>
#include <cctype>
#include <cmath>

#include "cprocess/param_db.hpp"

namespace cp {

namespace {

constexpr double NM = 1e-7;  // nm -> cm

std::string lower(std::string s) {
  std::transform(s.begin(), s.end(), s.begin(),
                 [](unsigned char c) { return std::tolower(c); });
  return s;
}

// Diffusivity prefactors/energies: Fair's vacancy model parameters for
// dopants in silicon (R.B. Fair, in "Impurity Doping Processes in Silicon",
// 1981; also the classic SUPREM defaults). Solid solubility fits and the
// implant range tables are coarse textbook-level approximations (LSS-style);
// the range table can always be overridden per-implant with rp=/drp=.
const std::vector<Dopant> kDopants = {
    {"boron", "B", DopType::acceptor, 5, 11.009,
     /*d0,e0*/ 0.037, 3.46, /*dm*/ 0, 0, /*dmm*/ 0, 0, /*dp*/ 0.72, 3.46,
     /*ss*/ 9.25e22, 0.73, /*fi*/ 1.0,
     {{10, 33 * NM, 17 * NM, -0.5, 3.5}, {20, 66 * NM, 28 * NM, -0.7, 4.0},
      {30, 99 * NM, 37 * NM, -0.8, 4.5}, {50, 161 * NM, 50 * NM, -1.0, 5.0},
      {80, 243 * NM, 63 * NM, -1.2, 6.0}, {100, 299 * NM, 71 * NM, -1.3, 6.5},
      {150, 420 * NM, 85 * NM, -1.4, 7.5}, {200, 531 * NM, 94 * NM, -1.5, 8.0}},
     /*seg_m0,seg_e*/ 6.0, 0.33, /*seg_h0,seg_he*/ 1.0e5, 2.0,
     /*dox0,eox*/ 1.23e-4, 3.39,
     /*dnit0,enit*/ 0.0, 0.0, /*dpoly0,epoly*/ 7.6, 3.46},
    {"phosphorus", "P", DopType::donor, 15, 30.974,
     3.85, 3.66, 4.44, 4.00, 44.2, 4.37, 0, 0,
     2.45e23, 0.62, /*fi*/ 1.0,
     {{10, 14 * NM, 7 * NM, -0.30, 3.2}, {20, 27 * NM, 13 * NM, -0.40, 3.4},
      {30, 42 * NM, 19 * NM, -0.45, 3.5}, {50, 68 * NM, 29 * NM, -0.55, 3.7},
      {80, 101 * NM, 40 * NM, -0.65, 4.0}, {100, 124 * NM, 45 * NM, -0.70, 4.2},
      {150, 190 * NM, 62 * NM, -0.80, 4.6}, {200, 254 * NM, 78 * NM, -0.90, 5.0}},
     /*seg_m0,seg_e*/ 10.0, 0.0, /*seg_h0,seg_he*/ 1.0e5, 2.0,
     /*dox0,eox*/ 0.19, 4.03,
     /*dnit0,enit*/ 0.0, 0.0, /*dpoly0,epoly*/ 40.0, 3.66},
    {"arsenic", "As", DopType::donor, 33, 74.922,
     0.066, 3.44, 12.0, 4.05, 0, 0, 0, 0,
     1.3e23, 0.66, /*fi*/ 0.4,
     {{10, 9 * NM, 4 * NM, -0.20, 3.1}, {20, 16 * NM, 7 * NM, -0.25, 3.2},
      {30, 23 * NM, 9 * NM, -0.30, 3.3}, {50, 34 * NM, 13 * NM, -0.35, 3.4},
      {80, 48 * NM, 18 * NM, -0.45, 3.6}, {100, 58 * NM, 21 * NM, -0.50, 3.7},
      {150, 85 * NM, 30 * NM, -0.60, 4.0}, {200, 110 * NM, 37 * NM, -0.70, 4.3}},
     /*seg_m0,seg_e*/ 10.0, 0.0, /*seg_h0,seg_he*/ 1.0e5, 2.0,
     /*dox0,eox*/ 3.7e-2, 3.70,
     /*dnit0,enit*/ 0.0, 0.0, /*dpoly0,epoly*/ 1.1, 3.44},
    {"antimony", "Sb", DopType::donor, 51, 120.90,
     0.214, 3.65, 15.0, 4.08, 0, 0, 0, 0,
     3.8e21, 0.56, /*fi*/ 0.1,
     {{10, 9 * NM, 3 * NM, -0.20, 3.1}, {30, 21 * NM, 7 * NM, -0.30, 3.3},
      {50, 31 * NM, 10 * NM, -0.35, 3.4}, {100, 53 * NM, 17 * NM, -0.50, 3.7},
      {200, 96 * NM, 29 * NM, -0.60, 4.0}},
     /*seg_m0,seg_e*/ 10.0, 0.0, /*seg_h0,seg_he*/ 1.0e5, 2.0,
     /*dox0,eox*/ 2.6e-2, 4.00,
     /*dnit0,enit*/ 0.0, 0.0, /*dpoly0,epoly*/ 5.3, 3.65},
};

}  // namespace

MatId material_id(const std::string& name) {
  const std::string q = lower(name);
  if (q == "silicon" || q == "si") return kMatSi;
  if (q == "oxide" || q == "sio2") return kMatOxide;
  if (q == "nitride" || q == "si3n4") return kMatNitride;
  if (q == "poly" || q == "polysilicon") return kMatPoly;
  return kMatGas;  // "gas" and anything unknown
}

const std::vector<Dopant>& dopant_table() { return kDopants; }

const Dopant* find_dopant(const std::string& name) {
  const std::string q = lower(name);
  for (const auto& d : kDopants)
    if (q == d.name || q == lower(d.symbol)) return &d;
  return nullptr;
}

double ni_si(double temp_k) {
  return 3.87e16 * std::pow(temp_k, 1.5) *
         std::exp(-0.605 / (kBoltzmannEv * temp_k));
}

double dopant_diffusivity(const Dopant& d, double temp_k, double n_over_ni) {
  const auto& P = ParamDB::instance();
  const double kt = kBoltzmannEv * temp_k;
  const double nni = std::max(n_over_ni, 1e-30);
  const double d0 = P.get(d.symbol + ".d0", d.d0);
  const double e0 = P.get(d.symbol + ".e0", d.e0);
  const double dm = P.get(d.symbol + ".dm", d.dm);
  const double em = P.get(d.symbol + ".em", d.em);
  const double dmm = P.get(d.symbol + ".dmm", d.dmm);
  const double emm = P.get(d.symbol + ".emm", d.emm);
  const double dp = P.get(d.symbol + ".dp", d.dp);
  const double ep = P.get(d.symbol + ".ep", d.ep);
  double dif = 0;
  if (d0 > 0) dif += d0 * std::exp(-e0 / kt);
  if (dm > 0) dif += dm * std::exp(-em / kt) * nni;
  if (dmm > 0) dif += dmm * std::exp(-emm / kt) * nni * nni;
  if (dp > 0) dif += dp * std::exp(-ep / kt) / nni;
  return dif;
}

double solid_solubility(const Dopant& d, double temp_k) {
  const auto& P = ParamDB::instance();
  const double ss_pre = P.get(d.symbol + ".ss_pre", d.ss_pre);
  const double ss_e = P.get(d.symbol + ".ss_e", d.ss_e);
  if (ss_pre <= 0) return 0;
  return ss_pre * std::exp(-ss_e / (kBoltzmannEv * temp_k));
}

double active_concentration(const Dopant& d, double conc, double temp_k) {
  const double css = solid_solubility(d, temp_k);
  return (css > 0) ? std::min(conc, css) : conc;
}

// Self-interstitial equilibrium concentration. Arrhenius fit giving ~1e13 cm^-3
// at 1000 C, rising toward ~1e15 near the melting point (cf. Bracht et al.).
double interstitial_cstar(double temp_k) {
  const auto& P = ParamDB::instance();
  const double pre = P.get("I.cstar_pre", 3.0e27);
  const double e = P.get("I.cstar_e", 3.7);
  return pre * std::exp(-e / (kBoltzmannEv * temp_k));
}

// Effective self-interstitial diffusivity (product D_I dominated by fast
// migration). ~1e-8 cm^2/s at 1000 C — the excess spreads and reaches the
// surface sink over seconds to minutes, setting the TED duration.
double interstitial_diffusivity(double temp_k) {
  const auto& P = ParamDB::instance();
  const double pre = P.get("I.d0", 5.0e-2);
  const double e = P.get("I.e0", 1.77);
  return pre * std::exp(-e / (kBoltzmannEv * temp_k));
}

// Bulk I-V recombination / trapping rate. The interstitial supersaturation
// relaxes toward equilibrium as exp(-k t) in the interior; combined with the
// surface sink this bounds the enhanced-diffusion transient. Calibrated so the
// transient lasts tens of seconds at typical anneal temperatures (tau ~ 1/k).
double interstitial_recomb_rate(double temp_k) {
  const auto& P = ParamDB::instance();
  const double pre = P.get("I.krec_pre", 2.0e4);
  const double e = P.get("I.krec_e", 1.4);
  return pre * std::exp(-e / (kBoltzmannEv * temp_k));
}

// Full point-defect model (P2-1). Vacancy defaults are chosen to be the same
// order of magnitude as the interstitial ones but distinct: C_V* uses a
// slightly lower activation energy (3.6 eV vs 3.7 eV) so vacancies are
// modestly more abundant at typical anneal temperatures, while D_V uses a
// higher prefactor-normalized but similar activation energy (1.8 eV vs
// 1.77 eV) and a 50x smaller prefactor, i.e. vacancies migrate distinctly
// slower than interstitials -- the physically-important ordering (I moves
// faster than V; both reach the same order of equilibrium concentration at
// high T) is what's calibrated here, not exact literature fits.
PointDefectParams point_defect_params(double temp_k, const ParamDB& db) {
  const double kt = kBoltzmannEv * temp_k;
  PointDefectParams p;
  const double ci_pre = db.get("pd.ci_star.pre", 3.0e27);
  const double ci_e = db.get("pd.ci_star.e", 3.7);
  p.ci_star = ci_pre * std::exp(-ci_e / kt);
  const double cv_pre = db.get("pd.cv_star.pre", 1.0e27);
  const double cv_e = db.get("pd.cv_star.e", 3.6);
  p.cv_star = cv_pre * std::exp(-cv_e / kt);
  const double di_pre = db.get("pd.di.pre", 5.0e-2);
  const double di_e = db.get("pd.di.e", 1.77);
  p.d_i = di_pre * std::exp(-di_e / kt);
  const double dv_pre = db.get("pd.dv.pre", 1.0e-3);
  const double dv_e = db.get("pd.dv.e", 1.8);
  p.d_v = dv_pre * std::exp(-dv_e / kt);

  const double a_si = 2.35e-8;  // cm, Si lattice capture radius
  const double kbulk_factor = db.get("pd.kbulk.factor", 1.0);
  p.k_bulk = kbulk_factor * 4.0 * M_PI * a_si * (p.d_i + p.d_v);
  const double ktrap_factor = db.get("pd.c311.ktrap.factor", 1.0);
  p.k_trap = ktrap_factor * 4.0 * M_PI * a_si * p.d_i;
  const double nu0 = db.get("pd.c311.nu0", 1e13);
  const double eb = db.get("pd.c311.eb", 3.6);
  p.k_emit = nu0 * std::exp(-eb / kt);
  return p;
}

double segregation_m(const Dopant& d, double temp_k) {
  const auto& P = ParamDB::instance();
  const double m0 = P.get(d.symbol + ".seg_m0", d.seg_m0);
  const double e = P.get(d.symbol + ".seg_e", d.seg_e);
  return m0 * std::exp(-e / (kBoltzmannEv * temp_k));
}

double segregation_h(const Dopant& d, double temp_k) {
  const auto& P = ParamDB::instance();
  const double h0 = P.get(d.symbol + ".seg_h0", d.seg_h0);
  const double he = P.get(d.symbol + ".seg_he", d.seg_he);
  return h0 * std::exp(-he / (kBoltzmannEv * temp_k));
}

double oxide_diffusivity(const Dopant& d, double temp_k) {
  const auto& P = ParamDB::instance();
  const double dox0 = P.get(d.symbol + ".dox0", d.dox0);
  const double eox = P.get(d.symbol + ".eox", d.eox);
  if (dox0 <= 0) return 0.0;
  return dox0 * std::exp(-eox / (kBoltzmannEv * temp_k));
}

double material_diffusivity(const Dopant& d, MatId mat, double temp_k,
                            double nni) {
  const auto& P = ParamDB::instance();
  switch (mat) {
    case kMatSi:
      return dopant_diffusivity(d, temp_k, nni);
    case kMatOxide:
      return oxide_diffusivity(d, temp_k);
    case kMatNitride: {
      const double dnit0 = P.get(d.symbol + ".dnit0", d.dnit0);
      const double enit = P.get(d.symbol + ".enit", d.enit);
      if (dnit0 <= 0) return 0.0;
      return dnit0 * std::exp(-enit / (kBoltzmannEv * temp_k));
    }
    case kMatPoly: {
      const double dpoly0 = P.get(d.symbol + ".dpoly0", d.dpoly0);
      const double epoly = P.get(d.symbol + ".epoly", d.epoly);
      if (dpoly0 <= 0) return 0.0;
      return dpoly0 * std::exp(-epoly / (kBoltzmannEv * temp_k));
    }
    case kMatGas:
    default:
      return 0.0;
  }
}

bool implant_moments(const Dopant& d, double energy_kev, double& rp,
                     double& drp, double& gamma, double& beta) {
  if (d.range.empty()) return false;
  const auto& t = d.range;
  if (energy_kev <= t.front()[0]) {
    rp = t.front()[1];
    drp = t.front()[2];
    gamma = t.front()[3];
    beta = t.front()[4];
    return true;
  }
  if (energy_kev >= t.back()[0]) {
    rp = t.back()[1];
    drp = t.back()[2];
    gamma = t.back()[3];
    beta = t.back()[4];
    return true;
  }
  for (std::size_t i = 1; i < t.size(); ++i) {
    if (energy_kev <= t[i][0]) {
      const double w = (std::log(energy_kev) - std::log(t[i - 1][0])) /
                       (std::log(t[i][0]) - std::log(t[i - 1][0]));
      rp = t[i - 1][1] + w * (t[i][1] - t[i - 1][1]);
      drp = t[i - 1][2] + w * (t[i][2] - t[i - 1][2]);
      gamma = t[i - 1][3] + w * (t[i][3] - t[i - 1][3]);
      beta = t[i - 1][4] + w * (t[i][4] - t[i - 1][4]);
      return true;
    }
  }
  return false;
}

bool implant_range(const Dopant& d, double energy_kev, double& rp, double& drp) {
  double gamma, beta;
  return implant_moments(d, energy_kev, rp, drp, gamma, beta);
}

}  // namespace cp
