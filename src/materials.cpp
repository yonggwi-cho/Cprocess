#include "cprocess/materials.hpp"

#include <algorithm>
#include <cctype>
#include <cmath>

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
     {{10, 33 * NM, 17 * NM}, {20, 66 * NM, 28 * NM}, {30, 99 * NM, 37 * NM},
      {50, 161 * NM, 50 * NM}, {80, 243 * NM, 63 * NM}, {100, 299 * NM, 71 * NM},
      {150, 420 * NM, 85 * NM}, {200, 531 * NM, 94 * NM}}},
    {"phosphorus", "P", DopType::donor, 15, 30.974,
     3.85, 3.66, 4.44, 4.00, 44.2, 4.37, 0, 0,
     2.45e23, 0.62, /*fi*/ 1.0,
     {{10, 14 * NM, 7 * NM}, {20, 27 * NM, 13 * NM}, {30, 42 * NM, 19 * NM},
      {50, 68 * NM, 29 * NM}, {80, 101 * NM, 40 * NM}, {100, 124 * NM, 45 * NM},
      {150, 190 * NM, 62 * NM}, {200, 254 * NM, 78 * NM}}},
    {"arsenic", "As", DopType::donor, 33, 74.922,
     0.066, 3.44, 12.0, 4.05, 0, 0, 0, 0,
     1.3e23, 0.66, /*fi*/ 0.4,
     {{10, 9 * NM, 4 * NM}, {20, 16 * NM, 7 * NM}, {30, 23 * NM, 9 * NM},
      {50, 34 * NM, 13 * NM}, {80, 48 * NM, 18 * NM}, {100, 58 * NM, 21 * NM},
      {150, 85 * NM, 30 * NM}, {200, 110 * NM, 37 * NM}}},
    {"antimony", "Sb", DopType::donor, 51, 120.90,
     0.214, 3.65, 15.0, 4.08, 0, 0, 0, 0,
     3.8e21, 0.56, /*fi*/ 0.1,
     {{10, 9 * NM, 3 * NM}, {30, 21 * NM, 7 * NM}, {50, 31 * NM, 10 * NM},
      {100, 53 * NM, 17 * NM}, {200, 96 * NM, 29 * NM}}},
};

}  // namespace

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
  const double kt = kBoltzmannEv * temp_k;
  const double nni = std::max(n_over_ni, 1e-30);
  double dif = 0;
  if (d.d0 > 0) dif += d.d0 * std::exp(-d.e0 / kt);
  if (d.dm > 0) dif += d.dm * std::exp(-d.em / kt) * nni;
  if (d.dmm > 0) dif += d.dmm * std::exp(-d.emm / kt) * nni * nni;
  if (d.dp > 0) dif += d.dp * std::exp(-d.ep / kt) / nni;
  return dif;
}

double solid_solubility(const Dopant& d, double temp_k) {
  if (d.ss_pre <= 0) return 0;
  return d.ss_pre * std::exp(-d.ss_e / (kBoltzmannEv * temp_k));
}

// Self-interstitial equilibrium concentration. Arrhenius fit giving ~1e13 cm^-3
// at 1000 C, rising toward ~1e15 near the melting point (cf. Bracht et al.).
double interstitial_cstar(double temp_k) {
  return 3.0e27 * std::exp(-3.7 / (kBoltzmannEv * temp_k));
}

// Effective self-interstitial diffusivity (product D_I dominated by fast
// migration). ~1e-8 cm^2/s at 1000 C — the excess spreads and reaches the
// surface sink over seconds to minutes, setting the TED duration.
double interstitial_diffusivity(double temp_k) {
  return 5.0e-2 * std::exp(-1.77 / (kBoltzmannEv * temp_k));
}

// Bulk I-V recombination / trapping rate. The interstitial supersaturation
// relaxes toward equilibrium as exp(-k t) in the interior; combined with the
// surface sink this bounds the enhanced-diffusion transient. Calibrated so the
// transient lasts tens of seconds at typical anneal temperatures (tau ~ 1/k).
double interstitial_recomb_rate(double temp_k) {
  return 2.0e4 * std::exp(-1.4 / (kBoltzmannEv * temp_k));
}

bool implant_range(const Dopant& d, double energy_kev, double& rp, double& drp) {
  if (d.range.empty()) return false;
  const auto& t = d.range;
  if (energy_kev <= t.front()[0]) {
    rp = t.front()[1];
    drp = t.front()[2];
    return true;
  }
  if (energy_kev >= t.back()[0]) {
    rp = t.back()[1];
    drp = t.back()[2];
    return true;
  }
  for (std::size_t i = 1; i < t.size(); ++i) {
    if (energy_kev <= t[i][0]) {
      const double w = (std::log(energy_kev) - std::log(t[i - 1][0])) /
                       (std::log(t[i][0]) - std::log(t[i - 1][0]));
      rp = t[i - 1][1] + w * (t[i][1] - t[i - 1][1]);
      drp = t[i - 1][2] + w * (t[i][2] - t[i - 1][2]);
      return true;
    }
  }
  return false;
}

}  // namespace cp
