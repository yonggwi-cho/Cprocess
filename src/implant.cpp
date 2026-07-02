#include "cprocess/implant.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace cp {

namespace {

struct Pearson4Table {
  double d_lo = 0, dz = 0;        // depth grid start / spacing [cm]
  std::vector<double> f;          // n points, integral(f * dz) = 1 (cm^-1)
  double eval(double depth) const {
    if (f.empty()) return 0.0;
    const double d_hi = d_lo + dz * (static_cast<double>(f.size()) - 1);
    if (depth < d_lo || depth > d_hi) return 0.0;
    const double t = (depth - d_lo) / dz;
    std::size_t i0 = static_cast<std::size_t>(t);
    if (i0 >= f.size() - 1) i0 = f.size() - 2;
    const double w = t - static_cast<double>(i0);
    return f[i0] + w * (f[i0 + 1] - f[i0]);
  }
};

// Returns false when the moments are invalid (caller falls back to Gauss).
bool build_pearson4(double rp, double drp, double gamma, double beta,
                    Pearson4Table& out) {
  const double A = 10.0 * beta - 12.0 * gamma * gamma - 18.0;
  if (A == 0.0) return false;
  const double b0 = -drp * drp * (4.0 * beta - 3.0 * gamma * gamma) / A;
  const double b1 = -gamma * drp * (beta + 3.0) / A;
  const double b2 = -(2.0 * beta - 3.0 * gamma * gamma - 6.0) / A;

  if (!(beta > 1.0 + gamma * gamma)) return false;
  if (!(b1 * b1 - 4.0 * b0 * b2 < 0.0)) return false;

  const double d_lo = std::max(0.0, rp - 6.0 * drp);
  const double d_hi = rp + 6.0 * drp;
  const int n = 2000;
  const double dz = (d_hi - d_lo) / (n - 1);

  // Pearson ODE d(ln f)/dz' = (z' - b1)/(b0 + b1 z' + b2 z'^2). With b0 < 0
  // the Gaussian limit (gamma=0, beta=3: b0=-dRp^2, b1=b2=0) gives
  // ln f = -z'^2/(2 dRp^2); integrating x^n times the ODE shows mean = 0
  // (i.e. Rp) and variance = -b0/(1+3 b2) = dRp^2 exactly.
  auto g = [&](double zp) {
    return (zp - b1) / (b0 + b1 * zp + b2 * zp * zp);
  };

  // Find grid index closest to Rp.
  int i_rp = static_cast<int>(std::round((rp - d_lo) / dz));
  if (i_rp < 0) i_rp = 0;
  if (i_rp >= n) i_rp = n - 1;

  std::vector<double> lnf(n, 0.0);
  lnf[i_rp] = 0.0;
  for (int i = i_rp + 1; i < n; ++i) {
    const double z_prev = (d_lo + (i - 1) * dz) - rp;
    const double z_cur = (d_lo + i * dz) - rp;
    lnf[i] = lnf[i - 1] + 0.5 * dz * (g(z_prev) + g(z_cur));
  }
  for (int i = i_rp - 1; i >= 0; --i) {
    const double z_next = (d_lo + (i + 1) * dz) - rp;
    const double z_cur = (d_lo + i * dz) - rp;
    lnf[i] = lnf[i + 1] - 0.5 * dz * (g(z_next) + g(z_cur));
  }

  std::vector<double> f(n);
  for (int i = 0; i < n; ++i) f[i] = std::exp(lnf[i]);

  // Trapezoidal-rule normalization so the grid integral of f is exactly 1.
  double s = 0.0;
  for (int i = 0; i < n - 1; ++i) s += 0.5 * dz * (f[i] + f[i + 1]);
  if (!(s > 0.0) || !std::isfinite(s)) return false;
  for (int i = 0; i < n; ++i) f[i] /= s;

  out.d_lo = d_lo;
  out.dz = dz;
  out.f = std::move(f);
  return true;
}

}  // namespace

double apply_implant(const Mesh& mesh, const std::vector<char>& mask,
                     const ImplantParams& p, std::vector<double>& conc) {
  if (!p.dopant) throw std::runtime_error("implant: no dopant");
  if (p.dose <= 0) throw std::runtime_error("implant: dose must be > 0");
  if (p.rp <= 0 || p.drp <= 0)
    throw std::runtime_error("implant: Rp and dRp must be > 0");

  const double drl = (p.drl > 0) ? p.drl : 0.8 * p.drp;
  const double ztop = mesh.bbox().hi.z;
  const double peak = p.dose / (std::sqrt(2.0 * M_PI) * p.drp);
  const double s2v = std::sqrt(2.0) * p.drp;
  const double s2l = std::sqrt(2.0) * drl;

  Pearson4Table tbl;
  const bool use_pearson = (p.profile == ImplantParams::Profile::pearson4) &&
                           build_pearson4(p.rp, p.drp, p.gamma, p.beta, tbl);

  conc.resize(mesh.cells.size(), 0.0);
  double atoms = 0;
  for (std::size_t ci = 0; ci < mesh.cells.size(); ++ci) {
    if (!mask.empty() && !mask[ci]) continue;
    const Vec3& c = mesh.cell_cent[ci];
    const double d = ztop - c.z;
    double v;
    if (use_pearson) {
      v = p.dose * tbl.eval(d);
    } else {
      v = peak * std::exp(-((d - p.rp) * (d - p.rp)) / (s2v * s2v));
    }
    if (p.has_window) {
      v *= 0.5 * (std::erf((c.x - p.x1) / s2l) - std::erf((c.x - p.x2) / s2l));
      v *= 0.5 * (std::erf((c.y - p.y1) / s2l) - std::erf((c.y - p.y2) / s2l));
    }
    if (v <= 0) continue;
    conc[ci] += v;
    atoms += v * mesh.cell_vol[ci];
  }
  return atoms;
}

}  // namespace cp
