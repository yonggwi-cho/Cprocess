#include "cprocess/mc_implant.hpp"

#include <array>
#include <cmath>
#include <cstdint>
#include <random>
#include <stdexcept>
#include <vector>

#ifdef _OPENMP
#  include <omp.h>
#else
static inline int omp_get_thread_num()  { return 0; }
static inline int omp_get_max_threads() { return 1; }
static inline void omp_set_num_threads(int) {}
#endif

#include "cprocess/cell_locator.hpp"

namespace cp {

namespace {

// Amorphous silicon target constants.
constexpr int    kZt      = 14;
constexpr double kMt      = 28.086;   // amu
constexpr double kNt      = 4.99e22;  // atoms/cm^3
constexpr double kEstopEv = 5.0;      // rest threshold [eV]
constexpr long long kChunk = 2048;    // ions per RNG stream / work unit

// Crystal silicon parameters (diamond cubic, a=5.431 Å).
constexpr double kAlatt  = 5.431e-8;  // lattice constant [cm]
constexpr double kU1     = 0.075e-8;  // 1-D thermal vibration amplitude [cm]
constexpr double kCL2    = 3.0;       // C_L^2 = 3 (Andersen-Feldman constant)
constexpr double ke2evcm = 14.3996e-8;// e^2 [eV·cm]

std::uint64_t splitmix64(std::uint64_t x) {
  x += 0x9E3779B97F4A7C15ull;
  x = (x ^ (x >> 30)) * 0xBF58476D1CE4E5B9ull;
  x = (x ^ (x >> 27)) * 0x94D049BB133111EBull;
  return x ^ (x >> 31);
}

struct Rng {
  std::mt19937_64 g;
  explicit Rng(std::uint64_t s) : g(s) {}
  double u() { return (g() >> 11) * (1.0 / 9007199254740992.0); }
};

// ZBL universal screening function and its derivative.
void zbl_screen(double x, double& phi, double& dphi) {
  static const double kC[4] = {0.18175, 0.50986, 0.28022, 0.02817};
  static const double kD[4] = {-3.19980, -0.94229, -0.40290, -0.20162};
  phi = 0; dphi = 0;
  for (int i = 0; i < 4; ++i) {
    const double e = kC[i] * std::exp(kD[i] * x);
    phi += e; dphi += kD[i] * e;
  }
}

// Reduced distance of closest approach: root of F(R)=1-V(R)/eps-(b/R)^2.
double closest_approach(double eps, double b) {
  double hi = std::max(b, 1.0);
  auto f = [&](double r) {
    double phi, dphi; zbl_screen(r, phi, dphi);
    return 1.0 - phi / (eps * r) - (b * b) / (r * r);
  };
  while (f(hi) < 0) hi *= 2.0;
  double lo = hi * 0.5;
  while (lo > 1e-12 && f(lo) > 0) { hi = lo; lo *= 0.5; }
  double r = 0.5 * (lo + hi);
  for (int it = 0; it < 60; ++it) {
    double phi, dphi; zbl_screen(r, phi, dphi);
    const double v = phi / r;
    const double fr = 1.0 - v / eps - (b * b) / (r * r);
    if (fr > 0) hi = r; else lo = r;
    const double dv = (dphi - v) / r;
    const double dfr = -dv / eps + 2.0 * b * b / (r * r * r);
    double rn = r - fr / dfr;
    if (!(rn > lo && rn < hi)) rn = 0.5 * (lo + hi);
    if (std::fabs(rn - r) < 1e-10 * r) return rn;
    r = rn;
  }
  return r;
}

// Exact classical scattering integral for the ZBL potential:
//   sin^2(theta/2) = cos^2(I),  I = int_0^1 g(w) dw,
// with R = R0/(1-w^2) substitution removing the sqrt singularity at R0.
// The w->0 limit is handled analytically via F'(R0) to avoid cancellation.
double integral_sin2(double eps, double b) {
  if (b <= 0) return 1.0;
  const double r0 = closest_approach(eps, b);
  double phi0, dphi0; zbl_screen(r0, phi0, dphi0);
  const double fp = (phi0 - r0 * dphi0) / (eps * r0 * r0) +
                    2.0 * b * b / (r0 * r0 * r0);
  const double lim0 = (fp > 0) ? 2.0*b / (r0 * std::sqrt(fp * r0)) : 0.0;

  const int n = 256;
  double sum = lim0 + 2.0 * (b / r0);  // g(0) and g(1)
  for (int i = 1; i < n; ++i) {
    const double w = static_cast<double>(i) / n;
    double g;
    if (w < 1e-3) {
      g = lim0;
    } else {
      const double r = r0 / (1.0 - w * w);
      double ph, dph; zbl_screen(r, ph, dph);
      const double fv = 1.0 - ph/(eps*r) - (b*b)/(r*r);
      g = 2.0 * w * (b/r0) / std::sqrt(std::max(fv, 1e-300));
    }
    sum += g * (i % 2 ? 4.0 : 2.0);
  }
  const double ci = std::cos(sum / (3.0 * n));
  return ci * ci;
}

// Tabulated sin^2(theta/2) on a log-log (eps, b) grid — corteo-style.
// Built once at startup; shared read-only across all threads.
struct ScatterTable {
  static constexpr int kNe = 160, kNb = 160;
  static constexpr double kEpsLo=1e-6, kEpsHi=1e3, kBLo=1e-4, kBHi=16.0;
  double le0, dle, lb0, dlb;
  std::vector<double> logs2;

  ScatterTable() : logs2(static_cast<std::size_t>(kNe) * kNb) {
    le0 = std::log(kEpsLo); dle = (std::log(kEpsHi) - le0) / (kNe - 1);
    lb0 = std::log(kBLo);   dlb = (std::log(kBHi)  - lb0) / (kNb - 1);
    for (int ie = 0; ie < kNe; ++ie) {
      const double eps = std::exp(le0 + ie * dle);
      for (int ib = 0; ib < kNb; ++ib) {
        const double b = std::exp(lb0 + ib * dlb);
        logs2[static_cast<std::size_t>(ie)*kNb + ib] =
            std::log(std::max(integral_sin2(eps, b), 1e-300));
      }
    }
  }

  double sample(double eps, double b) const {
    if (b <= kBLo) b = kBLo;
    if (b > kBHi)  b = kBHi;
    if (eps < kEpsLo) eps = kEpsLo;
    if (eps > kEpsHi) eps = kEpsHi;
    double u = (std::log(eps) - le0) / dle;
    double v = (std::log(b)   - lb0) / dlb;
    int i = std::min(static_cast<int>(u), kNe - 2);
    int j = std::min(static_cast<int>(v), kNb - 2);
    u -= i; v -= j;
    const double* t = logs2.data() + static_cast<std::size_t>(i)*kNb + j;
    const double l = (1-u)*((1-v)*t[0]+v*t[1]) + u*((1-v)*t[kNb]+v*t[kNb+1]);
    return std::exp(l);
  }
};

const ScatterTable& scatter_table() {
  static const ScatterTable table;
  return table;
}

}  // namespace

namespace mc {

double screening_length_cm(int z1, int z2) {
  return 0.46848e-8 / (std::pow(z1, 0.23) + std::pow(z2, 0.23));
}

double reduced_energy_per_ev(int z1, double m1, int z2, double m2) {
  return 0.032537 * m2 /
         (static_cast<double>(z1)*z2*(m1+m2)*
          (std::pow(z1, 0.23) + std::pow(z2, 0.23)));
}

double sin2_half_theta(double eps, double b) {
  return scatter_table().sample(eps, b);
}

}  // namespace mc

namespace {

Vec3 deflect(const Vec3& d, double cpsi, double spsi, double phi) {
  Vec3 e1 = (std::fabs(d.z) < 0.99) ? cross(d,{0,0,1}) : cross(d,{1,0,0});
  e1 = (1.0 / norm(e1)) * e1;
  const Vec3 e2 = cross(d, e1);
  Vec3 nd = cpsi*d + spsi*(std::cos(phi)*e1 + std::sin(phi)*e2);
  return (1.0 / norm(nd)) * nd;
}

// ---------------------------------------------------------------------------
// Crystal channeling model for Si (diamond cubic lattice)
// ---------------------------------------------------------------------------
// Lindhard-Robinson continuum string potential:
//   U(r) = (Z1 Z2 e²/d_row) × ln(1 + C_L² a_s²/r²)
// where d_row is the interatomic row spacing, a_s is the ZBL screening
// length, C_L²=3 (Andersen-Feldman), and U is evaluated in eV.
//
// Transverse energy: E_t = E sin²ψ (ignoring lateral oscillation energy).
// Channeling criterion: E_t < U_max = U(u₁), where u₁ = 0.075 Å is the
// 1-D thermal vibration amplitude.
//
// When channeling, the minimum impact parameter for nuclear collisions is
// raised from 0 to u₁, dramatically reducing nuclear stopping while
// keeping electronic stopping (slightly enhanced: ×1.2 for the higher
// electron density along an atomic row).
// ---------------------------------------------------------------------------

struct CrystalAxis {
  Vec3   dir;    // normalized crystal axis direction
  double u_max;  // Lindhard-Robinson U(u₁) [eV]; channel if E*sin²ψ < u_max
};

// Build channel axis list for a given ion species implanted into Si.
std::vector<CrystalAxis> build_crystal_axes(int z1, double a_s_cm) {
  auto u_at_u1 = [&](double d_row) {
    const double ratio = std::sqrt(kCL2) * a_s_cm / kU1;
    return ke2evcm * static_cast<double>(z1) * kZt / d_row *
           std::log(1.0 + ratio * ratio);
  };

  // Row spacings for Si diamond cubic:
  //   <100>: d = a  (atoms separated by one unit cell along [100])
  //   <110>: d = a/√2  (face diagonal row spacing)
  //   <111>: d = a√3/4  (body diagonal, closest-packed in diamond cubic)
  const double d100 = kAlatt;
  const double d110 = kAlatt * 0.70710678;  // / sqrt(2)
  const double d111 = kAlatt * 0.43301270;  // * sqrt(3)/4

  const double u100 = u_at_u1(d100);
  const double u110 = u_at_u1(d110);
  const double u111 = u_at_u1(d111);

  // 1/√2 ≈ 0.70711,  1/√3 ≈ 0.57735
  return {
    // <100> family (3 unique axes; ±axis gives same sin²ψ = 1-cos²ψ)
    {{1,0,0},u100}, {{0,1,0},u100}, {{0,0,1},u100},
    // <110> family (6 unique axes out of 12 ± directions)
    {{.70711,.70711,0},u110}, {{.70711,0,.70711},u110}, {{0,.70711,.70711},u110},
    {{.70711,-.70711,0},u110}, {{.70711,0,-.70711},u110}, {{0,.70711,-.70711},u110},
    // <111> family (4 unique axes)
    {{.57735,.57735,.57735},u111},  {{.57735,.57735,-.57735},u111},
    {{.57735,-.57735,.57735},u111}, {{-.57735,.57735,.57735},u111},
  };
}

// ---------------------------------------------------------------------------
// Walk parameters and ion transport
// ---------------------------------------------------------------------------

struct WalkParams {
  double e0_ev = 0;
  Vec3 dir0;
  double flight = 0;     // free-flight length [cm]
  double pmax = 0;       // max impact parameter [cm]
  double inv_a = 0;      // 1/screening length
  double eps_per_ev = 0;
  double tmax_fac = 0;
  double mass_ratio = 0;
  double els_fac = 0;    // electronic loss coefficient: dE = els_fac*sqrt(E)
  Vec3 lo, hi;
  bool wrap = false;
  double sx1, sx2, sy1, sy2;
  // Channeling (populated only when p.channeling == true)
  bool channeling = false;
  double u1_sq = 0;                     // u₁² [cm²]
  std::vector<CrystalAxis> chan_axes;   // crystal axes with precomputed U_max
};

enum class Fate { deposited, backscattered, transmitted, out_of_domain };

Fate walk_ion(const WalkParams& w, Rng& rng, Vec3& rest) {
  Vec3 pos = {w.sx1 + rng.u()*(w.sx2-w.sx1),
              w.sy1 + rng.u()*(w.sy2-w.sy1), w.hi.z};
  Vec3 dir = w.dir0;
  double e = w.e0_ev;

  auto wrap1 = [](double v, double lo, double hi) {
    const double span = hi - lo;
    v = std::fmod(v - lo, span);
    return lo + (v < 0 ? v + span : v);
  };

  while (true) {
    pos += w.flight * dir;
    if (pos.z > w.hi.z) return Fate::backscattered;
    if (pos.z < w.lo.z) return Fate::transmitted;
    if (pos.x < w.lo.x || pos.x > w.hi.x ||
        pos.y < w.lo.y || pos.y > w.hi.y) {
      if (!w.wrap) return Fate::out_of_domain;
      pos.x = wrap1(pos.x, w.lo.x, w.hi.x);
      pos.y = wrap1(pos.y, w.lo.y, w.hi.y);
    }

    // Channeling check: find nearest channel axis; test transverse energy.
    double eff_els = w.els_fac;
    double b_min_sq = 0.0;  // minimum impact parameter^2
    if (w.channeling) {
      double min_s2 = 1.0, u_max_best = 0.0;
      for (const auto& ax : w.chan_axes) {
        const double c = dot(dir, ax.dir);
        const double s2 = 1.0 - c * c;  // sin²ψ; symmetric in ±axis
        if (s2 < min_s2) { min_s2 = s2; u_max_best = ax.u_max; }
      }
      if (e * min_s2 < u_max_best) {
        // Ion is channeling along this axis.
        eff_els = w.els_fac * 1.2;  // enhanced electronic stopping in channel
        b_min_sq = w.u1_sq;         // nuclear collisions only beyond u₁
      }
    }

    // Continuous electronic loss over the free-flight segment.
    e -= eff_els * std::sqrt(e);
    if (e <= kEstopEv) break;

    // Binary nuclear collision: draw b uniformly in [b_min, pmax] area.
    const double p_sq = b_min_sq + rng.u() * (w.pmax*w.pmax - b_min_sq);
    const double p    = std::sqrt(p_sq);
    const double s2   = mc::sin2_half_theta(e * w.eps_per_ev, p * w.inv_a);
    e -= w.tmax_fac * e * s2;

    const double cth = 1.0 - 2.0*s2;
    const double sth = 2.0 * std::sqrt(std::max(0.0, s2*(1.0-s2)));
    const double psi = std::atan2(sth, cth + w.mass_ratio);
    dir = deflect(dir, std::cos(psi), std::sin(psi), 2.0*M_PI*rng.u());

    if (e <= kEstopEv) break;
  }
  rest = pos;
  return Fate::deposited;
}

}  // namespace

McImplantStats apply_mc_implant(const Mesh& mesh,
                                const std::vector<char>& silicon_mask,
                                const McImplantParams& p,
                                std::vector<double>& conc) {
  if (!p.dopant) throw std::runtime_error("mc implant: no dopant");
  if (p.dopant->z <= 0 || p.dopant->m <= 0)
    throw std::runtime_error("mc implant: species lacks Z/M data");
  if (p.dose <= 0) throw std::runtime_error("mc implant: dose must be > 0");
  if (p.energy_kev <= 0) throw std::runtime_error("mc implant: energy must be > 0");
  if (p.ions < 1) throw std::runtime_error("mc implant: ions must be >= 1");

  const int nc = static_cast<int>(mesh.cells.size());
  conc.resize(nc, 0.0);
  const CellLocator locator(mesh);
  const BBox bb = mesh.bbox();

  const int    z1 = p.dopant->z;
  const double m1 = p.dopant->m;

  WalkParams w;
  w.e0_ev     = p.energy_kev * 1000.0;
  w.flight    = std::pow(kNt, -1.0/3.0);
  w.pmax      = w.flight / std::sqrt(M_PI);
  w.inv_a     = 1.0 / mc::screening_length_cm(z1, kZt);
  w.eps_per_ev = mc::reduced_energy_per_ev(z1, m1, kZt, kMt);
  w.tmax_fac  = 4.0 * m1 * kMt / ((m1+kMt)*(m1+kMt));
  w.mass_ratio = m1 / kMt;
  // Lindhard-Scharff electronic stopping coefficient
  const double kls = 3.83e-15 * std::pow(z1, 7.0/6.0) * kZt /
                     (std::pow(std::pow(z1,2.0/3.0)+std::pow(kZt,2.0/3.0), 1.5) *
                      std::sqrt(m1));
  w.els_fac   = kNt * w.flight * kls / std::sqrt(1000.0);
  w.lo = bb.lo; w.hi = bb.hi;
  w.wrap = !p.has_window;
  w.sx1 = p.has_window ? p.x1 : bb.lo.x;
  w.sx2 = p.has_window ? p.x2 : bb.hi.x;
  w.sy1 = p.has_window ? p.y1 : bb.lo.y;
  w.sy2 = p.has_window ? p.y2 : bb.hi.y;
  const double tilt = p.tilt_deg * M_PI / 180.0;
  const double rot  = p.rotation_deg * M_PI / 180.0;
  w.dir0 = {std::sin(tilt)*std::cos(rot), std::sin(tilt)*std::sin(rot), -std::cos(tilt)};

  // Channeling setup: precompute crystal axes and U_max values.
  if (p.channeling) {
    w.channeling = true;
    w.u1_sq      = kU1 * kU1;
    w.chan_axes  = build_crystal_axes(z1, mc::screening_length_cm(z1, kZt));
  }

  scatter_table();  // build shared table before spawning workers

  const long long nchunks = (p.ions + kChunk - 1) / kChunk;

  // Determine thread count and configure OpenMP.
  int nthreads;
  if (p.threads > 0) {
    nthreads = p.threads;
    omp_set_num_threads(nthreads);
  } else {
    nthreads = omp_get_max_threads();
  }
  nthreads = std::max(1, std::min(nthreads, static_cast<int>(nchunks)));

  struct ChunkStat {
    long long dep=0, back=0, trans=0, ood=0, mask=0, unb=0;
    double dsum=0, d2sum=0;
  };
  std::vector<ChunkStat> cs(nchunks);
  std::vector<std::vector<std::uint32_t>> hits(
      nthreads, std::vector<std::uint32_t>(nc, 0));

  // Each chunk k gets a deterministic RNG seed derived from p.seed and k,
  // so results are bit-identical for any number of threads.
#pragma omp parallel for schedule(dynamic) num_threads(nthreads)
  for (long long k = 0; k < nchunks; ++k) {
    const int tid = omp_get_thread_num();
    std::vector<std::uint32_t>& h = hits[tid];
    Rng rng(splitmix64(p.seed ^ (0x9E3779B97F4A7C15ull * (k + 1))));
    ChunkStat& s = cs[k];
    const long long n = std::min<long long>(kChunk, p.ions - k * kChunk);
    for (long long i = 0; i < n; ++i) {
      Vec3 rest;
      switch (walk_ion(w, rng, rest)) {
        case Fate::backscattered: ++s.back; continue;
        case Fate::transmitted:   ++s.trans; continue;
        case Fate::out_of_domain: ++s.ood;  continue;
        case Fate::deposited:     break;
      }
      const int ci = locator.locate(rest);
      if (ci < 0) { ++s.unb; continue; }
      if (!silicon_mask.empty() && !silicon_mask[ci]) { ++s.mask; continue; }
      ++h[ci];
      ++s.dep;
      const double depth = bb.hi.z - rest.z;
      s.dsum  += depth;
      s.d2sum += depth * depth;
    }
  }

  // Sequential deterministic reduction (order of cs[] fixes floating-point).
  McImplantStats st;
  double dsum = 0, d2sum = 0;
  for (const auto& s : cs) {
    st.deposited     += s.dep;
    st.backscattered += s.back;
    st.transmitted   += s.trans;
    st.out_of_domain += s.ood;
    st.in_mask       += s.mask;
    st.unbinned      += s.unb;
    dsum  += s.dsum;
    d2sum += s.d2sum;
  }
  if (st.deposited > 0) {
    st.rp  = dsum / st.deposited;
    st.drp = std::sqrt(std::max(0.0, d2sum/st.deposited - st.rp*st.rp));
  }

  const double area = (w.sx2-w.sx1) * (w.sy2-w.sy1);
  st.weight = p.dose * area / static_cast<double>(p.ions);
  for (int i = 0; i < nc; ++i) {
    std::uint64_t total = 0;
    for (int t = 0; t < nthreads; ++t) total += hits[t][i];
    if (total > 0) conc[i] += total * st.weight / mesh.cell_vol[i];
  }
  return st;
}

}  // namespace cp
