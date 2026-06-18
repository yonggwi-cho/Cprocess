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
constexpr double kMt      = 28.086;    // amu
constexpr double kNt      = 4.99e22;   // atoms/cm^3
constexpr double kEstopEv = 5.0;       // ion rest threshold [eV]
constexpr long long kChunk = 2048;     // ions per RNG stream / work unit

// Crystal silicon (diamond cubic, a=5.431 Å).
constexpr double kAlatt  = 5.431e-8;   // lattice constant [cm]
constexpr double kU1     = 0.075e-8;   // 1-D thermal vibration amplitude [cm]
constexpr double kCL2    = 3.0;        // C_L^2=3 (Andersen-Feldman constant)
constexpr double ke2evcm = 14.3996e-8; // e^2 [eV·cm]

// Radiation damage.
constexpr double kEdSi     = 15.0;     // Si displacement threshold [eV]
constexpr double kNamorph  = 6.25e21;  // amorphization density [cm^-3] (~0.125*kNt)

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

// ZBL universal screening function and derivative.
void zbl_screen(double x, double& phi, double& dphi) {
  static const double kC[4] = {0.18175, 0.50986, 0.28022, 0.02817};
  static const double kD[4] = {-3.19980, -0.94229, -0.40290, -0.20162};
  phi = 0; dphi = 0;
  for (int i = 0; i < 4; ++i) {
    const double e = kC[i] * std::exp(kD[i] * x);
    phi += e; dphi += kD[i] * e;
  }
}

// Reduced distance of closest approach.
double closest_approach(double eps, double b) {
  double hi = std::max(b, 1.0);
  auto f = [&](double r) {
    double phi, dphi; zbl_screen(r, phi, dphi);
    return 1.0 - phi/(eps*r) - (b*b)/(r*r);
  };
  while (f(hi) < 0) hi *= 2.0;
  double lo = hi * 0.5;
  while (lo > 1e-12 && f(lo) > 0) { hi = lo; lo *= 0.5; }
  double r = 0.5*(lo+hi);
  for (int it = 0; it < 60; ++it) {
    double phi, dphi; zbl_screen(r, phi, dphi);
    const double v   = phi/r;
    const double fr  = 1.0 - v/eps - (b*b)/(r*r);
    if (fr > 0) hi = r; else lo = r;
    const double dv  = (dphi - v)/r;
    const double dfr = -dv/eps + 2.0*b*b/(r*r*r);
    double rn = r - fr/dfr;
    if (!(rn > lo && rn < hi)) rn = 0.5*(lo+hi);
    if (std::fabs(rn-r) < 1e-10*r) return rn;
    r = rn;
  }
  return r;
}

// Exact classical scattering integral for the ZBL potential;
// substitution R=R0/(1-w^2) removes the sqrt singularity at R0.
double integral_sin2(double eps, double b) {
  if (b <= 0) return 1.0;
  const double r0 = closest_approach(eps, b);
  double phi0, dphi0; zbl_screen(r0, phi0, dphi0);
  const double fp = (phi0 - r0*dphi0)/(eps*r0*r0) + 2.0*b*b/(r0*r0*r0);
  const double lim0 = (fp > 0) ? 2.0*b/(r0*std::sqrt(fp*r0)) : 0.0;
  const int n = 256;
  double sum = lim0 + 2.0*(b/r0);
  for (int i = 1; i < n; ++i) {
    const double w = static_cast<double>(i)/n;
    double g;
    if (w < 1e-3) {
      g = lim0;
    } else {
      const double r = r0/(1.0 - w*w);
      double ph, dph; zbl_screen(r, ph, dph);
      const double fv = 1.0 - ph/(eps*r) - (b*b)/(r*r);
      g = 2.0*w*(b/r0)/std::sqrt(std::max(fv, 1e-300));
    }
    sum += g * (i % 2 ? 4.0 : 2.0);
  }
  const double ci = std::cos(sum/(3.0*n));
  return ci*ci;
}

// Tabulated sin^2(theta/2) on a log-log (eps,b) grid — corteo-style.
// Immutable after construction; shared read-only across all worker threads.
struct ScatterTable {
  static constexpr int kNe=160, kNb=160;
  static constexpr double kEpsLo=1e-6, kEpsHi=1e3, kBLo=1e-4, kBHi=16.0;
  double le0, dle, lb0, dlb;
  std::vector<double> logs2;

  ScatterTable() : logs2(static_cast<std::size_t>(kNe)*kNb) {
    le0=std::log(kEpsLo); dle=(std::log(kEpsHi)-le0)/(kNe-1);
    lb0=std::log(kBLo);   dlb=(std::log(kBHi) -lb0)/(kNb-1);
    for (int ie=0; ie<kNe; ++ie) {
      const double eps=std::exp(le0+ie*dle);
      for (int ib=0; ib<kNb; ++ib) {
        const double b=std::exp(lb0+ib*dlb);
        logs2[static_cast<std::size_t>(ie)*kNb+ib]=
            std::log(std::max(integral_sin2(eps,b), 1e-300));
      }
    }
  }

  double sample(double eps, double b) const {
    if (b<=kBLo) b=kBLo; if (b>kBHi) b=kBHi;
    if (eps<kEpsLo) eps=kEpsLo; if (eps>kEpsHi) eps=kEpsHi;
    double u=(std::log(eps)-le0)/dle, v=(std::log(b)-lb0)/dlb;
    int i=std::min(static_cast<int>(u),kNe-2), j=std::min(static_cast<int>(v),kNb-2);
    u-=i; v-=j;
    const double* t=logs2.data()+static_cast<std::size_t>(i)*kNb+j;
    return std::exp((1-u)*((1-v)*t[0]+v*t[1])+u*((1-v)*t[kNb]+v*t[kNb+1]));
  }
};

const ScatterTable& scatter_table() { static const ScatterTable T; return T; }

}  // namespace

namespace mc {

double screening_length_cm(int z1, int z2) {
  return 0.46848e-8/(std::pow(z1,0.23)+std::pow(z2,0.23));
}
double reduced_energy_per_ev(int z1, double m1, int z2, double m2) {
  return 0.032537*m2/
         (static_cast<double>(z1)*z2*(m1+m2)*(std::pow(z1,0.23)+std::pow(z2,0.23)));
}
double sin2_half_theta(double eps, double b) { return scatter_table().sample(eps,b); }

}  // namespace mc

// Built-in stopping targets. Compounds use Bragg-rule effective Z/M and the
// mass-density-derived atomic number density N.
TargetMaterial target_silicon() {
  return {"Si", kZt, kMt, kNt, true};
}
TargetMaterial target_photoresist() {
  // DNQ-novolac ~ C-dominated organic, density ~1.2 g/cm^3. Effective carbon
  // atom: Z=6, M=12; N = rho/M * N_A = 1.2/12 * 6.022e23.
  return {"resist", 6, 12.0, 6.02e22, false};
}
TargetMaterial target_oxide() {
  // SiO2, density 2.2 g/cm^3, 3 atoms per 60 g/mol -> 20 g/mol per atom.
  // Effective Z = (14+8+8)/3 ~ 10, M ~ 20; N = 2.2/20 * 6.022e23.
  return {"SiO2", 10, 20.0, 6.62e22, false};
}
TargetMaterial target_vacuum() {
  // ~1000x rarefied air: large free flight, negligible stopping.
  return {"vacuum", 7, 14.0, 5.0e19, false};
}

namespace {

Vec3 deflect(const Vec3& d, double cpsi, double spsi, double phi) {
  Vec3 e1=(std::fabs(d.z)<0.99)?cross(d,{0,0,1}):cross(d,{1,0,0});
  e1=(1.0/norm(e1))*e1;
  const Vec3 e2=cross(d,e1);
  Vec3 nd=cpsi*d+spsi*(std::cos(phi)*e1+std::sin(phi)*e2);
  return (1.0/norm(nd))*nd;
}

// ---------------------------------------------------------------------------
// Crystal channeling: Lindhard-Robinson continuum string potential
//
//   U(r) = (Z1 Z2 e² / d_row) × ln(1 + C_L² a_s² / r²)
//
// C_L²=3 (Andersen-Feldman), a_s=ZBL screening length, d_row=row spacing.
// Channeling criterion per step: E·sin²ψ < U_max·(1-f_amor), where
// U_max=U(u₁) and f_amor=damage_density/kNamorph smoothly degrades the
// crystal as displacement damage accumulates.
// ---------------------------------------------------------------------------

struct CrystalAxis { Vec3 dir; double u_max; };

std::vector<CrystalAxis> build_crystal_axes(int z1, double a_s_cm) {
  auto u_at_u1 = [&](double d_row) {
    const double r = std::sqrt(kCL2) * a_s_cm / kU1;
    return ke2evcm * static_cast<double>(z1) * kZt / d_row * std::log(1.0+r*r);
  };
  const double d100 = kAlatt;
  const double d110 = kAlatt * 0.70710678;  // a/√2
  const double d111 = kAlatt * 0.43301270;  // a√3/4
  const double u100 = u_at_u1(d100), u110 = u_at_u1(d110), u111 = u_at_u1(d111);
  return {
    {{1,0,0},u100}, {{0,1,0},u100}, {{0,0,1},u100},
    {{.70711,.70711,0},u110}, {{.70711,0,.70711},u110}, {{0,.70711,.70711},u110},
    {{.70711,-.70711,0},u110}, {{.70711,0,-.70711},u110}, {{0,.70711,-.70711},u110},
    {{.57735,.57735,.57735},u111},  {{.57735,.57735,-.57735},u111},
    {{.57735,-.57735,.57735},u111}, {{-.57735,.57735,.57735},u111},
  };
}

// ---------------------------------------------------------------------------
// Shared damage accumulation
//
// All worker threads read and write this array via OpenMP atomic operations.
// No locking; relaxed visibility is acceptable: a slightly stale damage count
// only causes a minor over/under-estimate of the local amorphization fraction,
// which is a physically acceptable approximation.
// ---------------------------------------------------------------------------

struct DamageStore {
  std::vector<std::uint32_t> counts;  // displaced-atom counts per cell
  double weight;                      // real displaced atoms per simulated event
  const std::vector<double>* cell_vol;

  DamageStore() = default;
  DamageStore(int nc, double w, const std::vector<double>& cv)
      : counts(nc, 0), weight(w), cell_vol(&cv) {}

  // Add n displaced atoms to cell ci (called from parallel region).
  void add(int ci, std::uint32_t n) {
#pragma omp atomic
    counts[ci] += n;
  }

  // Amorphization fraction [0,1] at cell ci (read from parallel region).
  double f_amor(int ci) const {
    std::uint32_t cnt = 0;
#pragma omp atomic read
    cnt = counts[ci];
    const double dens = cnt * weight / (*cell_vol)[ci];
    return std::min(1.0, dens / kNamorph);
  }
};

// ---------------------------------------------------------------------------
// Walk parameters and ion transport
// ---------------------------------------------------------------------------

// Precomputed ion-in-material stopping constants for one (dopant, target) pair.
struct MatConstants {
  double flight = 0, pmax = 0, inv_a = 0;
  double eps_per_ev = 0, tmax_fac = 0, mass_ratio = 0, els_fac = 0;
  bool crystal_si = false;
  std::vector<CrystalAxis> chan_axes;  // populated only for crystal Si
};

MatConstants make_mat_constants(int z1, double m1, const TargetMaterial& t) {
  MatConstants c;
  const int z2 = t.z;
  const double m2 = t.m, N = t.n;
  c.flight     = std::pow(N, -1.0/3.0);
  c.pmax       = c.flight / std::sqrt(M_PI);
  c.inv_a      = 1.0 / mc::screening_length_cm(z1, z2);
  c.eps_per_ev = mc::reduced_energy_per_ev(z1, m1, z2, m2);
  c.tmax_fac   = 4.0*m1*m2/((m1+m2)*(m1+m2));
  c.mass_ratio = m1/m2;
  const double kls = 3.83e-15*std::pow(z1,7.0/6.0)*z2/
      (std::pow(std::pow(z1,2.0/3.0)+std::pow(z2,2.0/3.0),1.5)*std::sqrt(m1));
  c.els_fac    = N*c.flight*kls/std::sqrt(1000.0);
  c.crystal_si = t.crystal_si;
  return c;
}

struct WalkParams {
  double e0_ev = 0;
  Vec3 dir0;
  Vec3 lo, hi;
  bool wrap = false;
  double sx1, sx2, sy1, sy2;
  // Per-material stopping constants; mats[0] is the substrate.
  std::vector<MatConstants> mats;
  bool single_material = true;             // true -> use mats[0], skip per-step locate
  const std::vector<int>* cell_material = nullptr;
  int default_mat = 0;                      // material when a cell is not found
  // Channeling
  bool channeling = false;
  double u1_sq = 0;
  // Damage — non-null when channeling is enabled
  DamageStore* damage = nullptr;
  const CellLocator* locator = nullptr;
};

enum class Fate { deposited, backscattered, transmitted, out_of_domain };

Fate walk_ion(const WalkParams& w, Rng& rng, Vec3& rest) {
  Vec3 pos = {w.sx1+rng.u()*(w.sx2-w.sx1), w.sy1+rng.u()*(w.sy2-w.sy1), w.hi.z};
  Vec3 dir = w.dir0;
  double e = w.e0_ev;

  auto wrap1 = [](double v, double lo, double hi) {
    const double span=hi-lo; v=std::fmod(v-lo,span);
    return lo+(v<0?v+span:v);
  };

  while (true) {
    // Material the ion is currently traversing. In single-material mode this is
    // always the substrate (no locate, preserving speed and determinism).
    int mid = w.default_mat;
    if (!w.single_material) {
      const int ci = w.locator->locate(pos);
      if (ci >= 0 && w.cell_material) mid = (*w.cell_material)[ci];
    }
    const MatConstants& mc_ = w.mats[mid];

    pos += mc_.flight * dir;
    if (pos.z > w.hi.z) return Fate::backscattered;
    if (pos.z < w.lo.z) return Fate::transmitted;
    if (pos.x<w.lo.x||pos.x>w.hi.x||pos.y<w.lo.y||pos.y>w.hi.y) {
      if (!w.wrap) return Fate::out_of_domain;
      pos.x=wrap1(pos.x,w.lo.x,w.hi.x);
      pos.y=wrap1(pos.y,w.lo.y,w.hi.y);
    }

    // Channeling only inside crystalline silicon.
    double eff_els = mc_.els_fac;
    double b_min_sq = 0.0;
    if (w.channeling && mc_.crystal_si) {
      double min_s2 = 1.0, u_max_best = 0.0;
      for (const auto& ax : mc_.chan_axes) {
        const double c = dot(dir, ax.dir);
        const double s2 = 1.0 - c*c;
        if (s2 < min_s2) { min_s2=s2; u_max_best=ax.u_max; }
      }
      // Reduce effective U_max by amorphization fraction at current cell.
      double f = 0.0;
      if (w.damage) {
        const int ci = w.locator->locate(pos);
        if (ci >= 0) f = w.damage->f_amor(ci);
      }
      const double u_eff = u_max_best * (1.0 - f);
      if (e * min_s2 < u_eff) {
        eff_els = mc_.els_fac * 1.2;
        b_min_sq = w.u1_sq;
      }
    }

    // Continuous electronic loss over the free-flight segment.
    e -= eff_els * std::sqrt(e);
    if (e <= kEstopEv) break;

    // Binary nuclear collision.
    const double p_sq = b_min_sq + rng.u()*(mc_.pmax*mc_.pmax - b_min_sq);
    const double p    = std::sqrt(p_sq);
    const double s2   = mc::sin2_half_theta(e * mc_.eps_per_ev, p * mc_.inv_a);
    const double t_recoil = mc_.tmax_fac * e * s2;  // recoil energy [eV]
    e -= t_recoil;

    // Kinchin-Pease damage: accumulate displaced atoms (crystal Si only).
    if (w.damage && mc_.crystal_si && t_recoil >= kEdSi) {
      const int ci = w.locator->locate(pos);
      if (ci >= 0) {
        const std::uint32_t ndis = (t_recoil < 2.0*kEdSi)
            ? 1u
            : static_cast<std::uint32_t>(t_recoil / (2.0*kEdSi));
        w.damage->add(ci, ndis);
      }
    }

    const double cth = 1.0 - 2.0*s2;
    const double sth = 2.0*std::sqrt(std::max(0.0, s2*(1.0-s2)));
    const double psi = std::atan2(sth, cth+mc_.mass_ratio);
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
                                std::vector<double>& conc,
                                std::vector<double>* damage_conc) {
  if (!p.dopant) throw std::runtime_error("mc implant: no dopant");
  if (p.dopant->z<=0||p.dopant->m<=0)
    throw std::runtime_error("mc implant: species lacks Z/M data");
  if (p.dose<=0) throw std::runtime_error("mc implant: dose must be > 0");
  if (p.energy_kev<=0) throw std::runtime_error("mc implant: energy must be > 0");
  if (p.ions<1) throw std::runtime_error("mc implant: ions must be >= 1");

  const int nc = static_cast<int>(mesh.cells.size());
  conc.resize(nc, 0.0);
  const CellLocator locator(mesh);
  const BBox bb = mesh.bbox();

  const int    z1 = p.dopant->z;
  const double m1 = p.dopant->m;

  WalkParams w;
  w.e0_ev     = p.energy_kev * 1000.0;

  // Build per-material stopping constants. With no material table the domain is
  // a single crystalline-Si target (legacy behaviour).
  if (p.cell_material && !p.material_table.empty()) {
    w.single_material = false;
    w.cell_material   = p.cell_material;
    w.mats.reserve(p.material_table.size());
    for (const auto& t : p.material_table) w.mats.push_back(make_mat_constants(z1, m1, t));
  } else {
    w.single_material = true;
    w.mats.push_back(make_mat_constants(z1, m1, target_silicon()));
  }
  w.default_mat = 0;

  w.lo=bb.lo; w.hi=bb.hi;
  w.wrap = !p.has_window;
  w.sx1=p.has_window?p.x1:bb.lo.x; w.sx2=p.has_window?p.x2:bb.hi.x;
  w.sy1=p.has_window?p.y1:bb.lo.y; w.sy2=p.has_window?p.y2:bb.hi.y;
  const double tilt=p.tilt_deg*M_PI/180.0, rot=p.rotation_deg*M_PI/180.0;
  w.dir0={std::sin(tilt)*std::cos(rot), std::sin(tilt)*std::sin(rot), -std::cos(tilt)};

  // Precompute ion weight for the damage store.
  const double area = (w.sx2-w.sx1)*(w.sy2-w.sy1);
  const double weight = p.dose * area / static_cast<double>(p.ions);

  // Channeling setup. Crystal axes are attached to every crystalline-Si
  // material so the ion channels only while inside silicon.
  DamageStore dmg_store;
  if (p.channeling) {
    w.channeling = true;
    w.u1_sq      = kU1*kU1;
    for (auto& mc_ : w.mats)
      if (mc_.crystal_si)
        mc_.chan_axes = build_crystal_axes(z1, mc::screening_length_cm(z1, kZt));
    dmg_store    = DamageStore(nc, weight, mesh.cell_vol);
    w.damage     = &dmg_store;
    w.locator    = &locator;
  }
  // Multi-material transport needs the locator even without channeling.
  if (!w.single_material && !w.locator) w.locator = &locator;

  scatter_table();  // ensure table is built before spawning workers

  const long long nchunks = (p.ions + kChunk - 1) / kChunk;
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

  // Each chunk k uses a fixed seed derived from p.seed, ensuring results are
  // bit-identical for any number of threads (amorphous mode).
  // Channeling with live damage feedback is non-deterministic across thread
  // counts by design: scheduling affects the order in which damage accumulates.
#pragma omp parallel for schedule(dynamic) num_threads(nthreads)
  for (long long k = 0; k < nchunks; ++k) {
    const int tid = omp_get_thread_num();
    std::vector<std::uint32_t>& h = hits[tid];
    Rng rng(splitmix64(p.seed ^ (0x9E3779B97F4A7C15ull*(k+1))));
    ChunkStat& s = cs[k];
    const long long n = std::min<long long>(kChunk, p.ions - k*kChunk);
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
      s.d2sum += depth*depth;
    }
  }

  // Sequential reduction (fixed order preserves floating-point for dsum/d2sum).
  McImplantStats st;
  double dsum=0, d2sum=0;
  for (const auto& s : cs) {
    st.deposited     += s.dep;  st.backscattered += s.back;
    st.transmitted   += s.trans; st.out_of_domain += s.ood;
    st.in_mask       += s.mask;  st.unbinned      += s.unb;
    dsum += s.dsum; d2sum += s.d2sum;
  }
  if (st.deposited > 0) {
    st.rp  = dsum / st.deposited;
    st.drp = std::sqrt(std::max(0.0, d2sum/st.deposited - st.rp*st.rp));
  }

  st.weight = weight;
  for (int i = 0; i < nc; ++i) {
    std::uint64_t total = 0;
    for (int t = 0; t < nthreads; ++t) total += hits[t][i];
    if (total > 0) conc[i] += total * st.weight / mesh.cell_vol[i];
  }

  // Optionally expose damage profile (displaced-atom density [cm^-3]).
  if (damage_conc && p.channeling) {
    damage_conc->assign(nc, 0.0);
    for (int i = 0; i < nc; ++i)
      (*damage_conc)[i] = dmg_store.counts[i] * weight / mesh.cell_vol[i];
  }

  return st;
}

}  // namespace cp
