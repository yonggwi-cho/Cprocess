#include "cprocess/diffusion.hpp"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <limits>
#include <ostream>
#include <stdexcept>

namespace cp {

namespace {
const double kNaN = std::numeric_limits<double>::quiet_NaN();
}

DiffusionSolver::DiffusionSolver(const Mesh& mesh, std::vector<char> solve_mask,
                                 std::ostream* log)
    : mesh_(mesh), mask_(std::move(solve_mask)), log_(log) {
  if (mask_.empty()) mask_.assign(mesh_.cells.size(), 1);
  if (mask_.size() != mesh_.cells.size())
    throw std::runtime_error("diffusion: mask size mismatch");
  build();
}

void DiffusionSolver::build() {
  const int nc = static_cast<int>(mesh_.cells.size());
  const int nf = static_cast<int>(mesh_.faces.size());
  fg_.assign(nf, FGeom{});
  clamped_faces_ = 0;

  for (int fi = 0; fi < nf; ++fi) {
    const Face& f = mesh_.faces[fi];
    const bool mo = mask_[f.owner];
    const bool mn = f.neigh >= 0 && mask_[f.neigh];
    FGeom& g = fg_[fi];
    if (mo && mn) {
      g.kind = kInternal;
      const Vec3 d = mesh_.cell_cent[f.neigh] - mesh_.cell_cent[f.owner];
      const double sn = norm(f.S), dn = norm(d);
      double sd = dot(f.S, d);
      // Guard against severely non-orthogonal faces; the clamp keeps the
      // matrix an M-matrix at the cost of local accuracy.
      if (sd < 0.05 * sn * dn) {
        sd = 0.05 * sn * dn;
        ++clamped_faces_;
      }
      g.g = dot(f.S, f.S) / sd;
      g.k = f.S - g.g * d;
      const Vec3 nh = (1.0 / sn) * f.S;
      g.delP = std::max(std::fabs(dot(f.c - mesh_.cell_cent[f.owner], nh)), 1e-3 * dn);
      g.delN = std::max(std::fabs(dot(mesh_.cell_cent[f.neigh] - f.c, nh)), 1e-3 * dn);
      g.wP = g.delN / (g.delP + g.delN);
    } else if (mo) {
      g.kind = kBoundOwner;
      const Vec3 db = f.c - mesh_.cell_cent[f.owner];
      const double sn = norm(f.S), dn = norm(db);
      const double sd = std::max(dot(f.S, db), 0.05 * sn * dn + 1e-300);
      g.gb = dot(f.S, f.S) / sd;
    } else if (mn) {
      g.kind = kBoundNeigh;
    }
  }

  // CSR pattern: one row per cell; masked cells keep an identity row.
  std::vector<std::vector<int>> nb(nc);
  for (int fi = 0; fi < nf; ++fi) {
    if (fg_[fi].kind != kInternal) continue;
    const Face& f = mesh_.faces[fi];
    nb[f.owner].push_back(f.neigh);
    nb[f.neigh].push_back(f.owner);
  }
  A_.n = nc;
  A_.ptr.assign(nc + 1, 0);
  for (int i = 0; i < nc; ++i) {
    auto& v = nb[i];
    v.push_back(i);
    std::sort(v.begin(), v.end());
    v.erase(std::unique(v.begin(), v.end()), v.end());
    A_.ptr[i + 1] = A_.ptr[i] + static_cast<int>(v.size());
  }
  A_.col.resize(A_.ptr[nc]);
  for (int i = 0; i < nc; ++i)
    std::copy(nb[i].begin(), nb[i].end(), A_.col.begin() + A_.ptr[i]);
  A_.val.assign(A_.ptr[nc], 0.0);

  diag_.resize(nc);
  for (int i = 0; i < nc; ++i) diag_[i] = A_.find(i, i);
  fslot_.assign(nf, {-1, -1});
  for (int fi = 0; fi < nf; ++fi) {
    if (fg_[fi].kind != kInternal) continue;
    const Face& f = mesh_.faces[fi];
    fslot_[fi] = {A_.find(f.owner, f.neigh), A_.find(f.neigh, f.owner)};
  }

  if (log_ && clamped_faces_ > 0)
    *log_ << "[diffuse] warning: " << clamped_faces_
          << " faces with poor orthogonality were clamped\n";
}

namespace {

// Solves the SPD-ish 3x3 system G x = b with partial pivoting; G is
// Tikhonov-regularized by the caller, so a tiny pivot means a genuinely
// missing direction (the component stays ~0).
Vec3 solve3(double G[3][3], const Vec3& b) {
  double a[3][4] = {{G[0][0], G[0][1], G[0][2], b.x},
                    {G[1][0], G[1][1], G[1][2], b.y},
                    {G[2][0], G[2][1], G[2][2], b.z}};
  for (int k = 0; k < 3; ++k) {
    int piv = k;
    for (int i = k + 1; i < 3; ++i)
      if (std::fabs(a[i][k]) > std::fabs(a[piv][k])) piv = i;
    if (piv != k)
      for (int j = 0; j < 4; ++j) std::swap(a[k][j], a[piv][j]);
    if (a[k][k] == 0.0) return {0, 0, 0};
    for (int i = k + 1; i < 3; ++i) {
      const double m = a[i][k] / a[k][k];
      for (int j = k; j < 4; ++j) a[i][j] -= m * a[k][j];
    }
  }
  Vec3 x;
  x.z = a[2][3] / a[2][2];
  x.y = (a[1][3] - a[1][2] * x.z) / a[1][1];
  x.x = (a[0][3] - a[0][1] * x.y - a[0][2] * x.z) / a[0][0];
  return x;
}

}  // namespace

void DiffusionSolver::gradients(const std::vector<double>& c,
                                const std::vector<double>& bcface,
                                std::vector<Vec3>& grad) const {
  // Weighted least-squares gradient over face neighbours (exact for linear
  // fields on any mesh). Boundary closure: Dirichlet faces contribute a
  // ghost point at the face centroid; zero-flux faces contribute the
  // constraint grad.n = 0, which both regularizes boundary cells and
  // encodes the physical condition.
  const int nc = static_cast<int>(mesh_.cells.size());
  std::vector<std::array<double, 9>> G(nc, {0, 0, 0, 0, 0, 0, 0, 0, 0});
  grad.assign(nc, Vec3{});

  auto add_row = [&](int cell, Vec3 d, double dc) {
    const double w = 1.0 / std::max(dot(d, d), 1e-300);
    auto& g = G[cell];
    g[0] += w * d.x * d.x; g[1] += w * d.x * d.y; g[2] += w * d.x * d.z;
    g[4] += w * d.y * d.y; g[5] += w * d.y * d.z; g[8] += w * d.z * d.z;
    grad[cell] += (w * dc) * d;  // accumulate rhs in `grad`
  };

  for (std::size_t fi = 0; fi < mesh_.faces.size(); ++fi) {
    const Face& f = mesh_.faces[fi];
    switch (fg_[fi].kind) {
      case kInternal: {
        const Vec3 d = mesh_.cell_cent[f.neigh] - mesh_.cell_cent[f.owner];
        const double dc = c[f.neigh] - c[f.owner];
        add_row(f.owner, d, dc);
        add_row(f.neigh, -1.0 * d, -dc);
        break;
      }
      case kBoundOwner: {
        const Vec3 d = f.c - mesh_.cell_cent[f.owner];
        if (!std::isnan(bcface[fi])) {
          add_row(f.owner, d, bcface[fi] - c[f.owner]);
        } else {
          const Vec3 nh = (1.0 / norm(f.S)) * f.S;
          add_row(f.owner, dot(d, nh) * nh, 0.0);
        }
        break;
      }
      case kBoundNeigh: {
        const Vec3 d = f.c - mesh_.cell_cent[f.neigh];
        const Vec3 nh = (1.0 / norm(f.S)) * f.S;
        add_row(f.neigh, dot(d, nh) * nh, 0.0);
        break;
      }
      default:
        break;
    }
  }

  for (int i = 0; i < nc; ++i) {
    if (!mask_[i]) { grad[i] = Vec3{}; continue; }
    auto& g = G[i];
    const double tr = g[0] + g[4] + g[8];
    const double eps = 1e-10 * (tr > 0 ? tr : 1.0);
    double m[3][3] = {{g[0] + eps, g[1], g[2]},
                      {g[1], g[4] + eps, g[5]},
                      {g[2], g[5], g[8] + eps}};
    grad[i] = solve3(m, grad[i]);
  }
}

void DiffusionSolver::run(std::vector<SpeciesField>& fields,
                          const std::vector<DirichletBC>& bcs,
                          const DiffuseOpts& o) {
  const int nc = static_cast<int>(mesh_.cells.size());
  const int nf = static_cast<int>(mesh_.faces.size());
  const int ns = static_cast<int>(fields.size());
  if (ns == 0) return;
  for (auto& sf : fields) {
    if (!sf.dopant || !sf.conc) throw std::runtime_error("diffusion: bad field");
    sf.conc->resize(nc, 0.0);
  }
  std::sort(fields.begin(), fields.end(),
            [](const SpeciesField& a, const SpeciesField& b) {
              return a.dopant->symbol < b.dopant->symbol;
            });

  const double dt0 = (o.dt > 0) ? std::min(o.dt, o.time) : o.time / 50.0;
  const int nsteps = static_cast<int>(std::ceil(o.time / dt0 - 1e-12));
  const double ni = ni_si(o.temp);

  // Per-species Dirichlet value per boundary face (NaN = no condition).
  std::vector<std::vector<double>> bcface(ns, std::vector<double>(nf, kNaN));
  for (int s = 0; s < ns; ++s)
    for (const auto& bc : bcs) {
      if (bc.species != fields[s].dopant->symbol || bc.patch < 0) continue;
      for (int fi = 0; fi < nf; ++fi)
        if (mesh_.faces[fi].patch == bc.patch && fg_[fi].kind == kBoundOwner &&
            mesh_.faces[fi].neigh < 0)
          bcface[s][fi] = bc.conc;
    }

  std::vector<double> mass0(ns, 0.0);
  for (int s = 0; s < ns; ++s)
    for (int i = 0; i < nc; ++i)
      mass0[s] += (*fields[s].conc)[i] * mesh_.cell_vol[i];

  std::vector<std::vector<double>> cold(ns), dcell(ns, std::vector<double>(nc, 0.0));
  std::vector<double> nni(nc, 1.0), rhs(nc), x;
  std::vector<Vec3> grad;

  double t = 0;
  for (int step = 0; step < nsteps; ++step) {
    const double dt = std::min(dt0, o.time - t);
    for (int s = 0; s < ns; ++s) cold[s] = *fields[s].conc;

    double cmax = 0;
    for (int s = 0; s < ns; ++s)
      for (int i = 0; i < nc; ++i) cmax = std::max(cmax, (*fields[s].conc)[i]);
    const double cfloor = 1e-3 * std::max(cmax, 1.0);

    int picard = 0, lin_iters = 0;
    double maxrel = 0;
    for (picard = 1; picard <= o.max_picard; ++picard) {
      // n/ni from charge neutrality with the current iterate of all species.
      for (int i = 0; i < nc; ++i) {
        double nnet = 0;
        for (int s = 0; s < ns; ++s) {
          const double c = (*fields[s].conc)[i];
          nnet += (fields[s].dopant->type == DopType::donor) ? c : -c;
        }
        const double cc = nnet / (2.0 * ni);
        nni[i] = cc + std::sqrt(cc * cc + 1.0);
      }
      // Per-cell diffusivity, optionally with field enhancement.
      for (int s = 0; s < ns; ++s) {
        const Dopant& dp = *fields[s].dopant;
        for (int i = 0; i < nc; ++i) {
          if (!mask_[i]) { dcell[s][i] = 0; continue; }
          double dv = dopant_diffusivity(dp, o.temp, nni[i]);
          if (o.field_enh) {
            const bool ntype = nni[i] >= 1.0;
            if ((dp.type == DopType::donor && ntype) ||
                (dp.type == DopType::acceptor && !ntype)) {
              const double cc = 0.5 * (nni[i] - 1.0 / nni[i]);  // = Nnet/(2 ni)
              dv *= 1.0 + std::fabs(cc) / std::sqrt(cc * cc + 1.0);
            }
          }
          dcell[s][i] = dv;
        }
      }

      maxrel = 0;
      for (int s = 0; s < ns; ++s) {
        std::vector<double>& c = *fields[s].conc;
        gradients(c, bcface[s], grad);

        std::fill(A_.val.begin(), A_.val.end(), 0.0);
        for (int i = 0; i < nc; ++i) {
          if (!mask_[i]) {
            A_.val[diag_[i]] = 1.0;
            rhs[i] = c[i];
          } else {
            A_.val[diag_[i]] = mesh_.cell_vol[i] / dt;
            rhs[i] = mesh_.cell_vol[i] / dt * cold[s][i];
          }
        }
        for (int fi = 0; fi < nf; ++fi) {
          const FGeom& g = fg_[fi];
          const Face& f = mesh_.faces[fi];
          if (g.kind == kInternal) {
            const double dP = dcell[s][f.owner], dN = dcell[s][f.neigh];
            if (dP <= 0 || dN <= 0) continue;
            const double dh = (g.delP + g.delN) / (g.delP / dP + g.delN / dN);
            const double tf = dh * g.g;
            A_.val[diag_[f.owner]] += tf;
            A_.val[diag_[f.neigh]] += tf;
            A_.val[fslot_[fi][0]] -= tf;
            A_.val[fslot_[fi][1]] -= tf;
            if (o.nonortho) {
              const Vec3 gf = g.wP * grad[f.owner] + (1.0 - g.wP) * grad[f.neigh];
              const double corr = dh * dot(gf, g.k);
              rhs[f.owner] += corr;
              rhs[f.neigh] -= corr;
            }
          } else if (g.kind == kBoundOwner && !std::isnan(bcface[s][fi])) {
            const double tb = dcell[s][f.owner] * g.gb;
            A_.val[diag_[f.owner]] += tb;
            rhs[f.owner] += tb * bcface[s][fi];
          }
        }

        x = c;  // warm start
        SolveResult sr = cg_jacobi(A_, rhs, x, o.lin_rtol, o.lin_maxit);
        if (!sr.converged) {
          x = c;
          sr = bicgstab_jacobi(A_, rhs, x, o.lin_rtol, o.lin_maxit);
          if (!sr.converged)
            throw std::runtime_error(
                "diffusion: linear solver failed (resid=" +
                std::to_string(sr.resid) + ")");
        }
        lin_iters += sr.iters;
        for (int i = 0; i < nc; ++i)
          maxrel = std::max(maxrel,
                            std::fabs(x[i] - c[i]) / (std::fabs(x[i]) + cfloor));
        c = x;
      }
      if (maxrel < o.picard_tol) break;
    }

    // The deferred correction can produce small negative undershoots.
    for (int s = 0; s < ns; ++s)
      for (double& v : *fields[s].conc) v = std::max(v, 0.0);

    t += dt;
    if (o.verbosity >= 1 && log_) {
      const int every = std::max(1, nsteps / 10);
      if (step % every == 0 || step == nsteps - 1) {
        char buf[160];
        std::snprintf(buf, sizeof(buf),
                      "[diffuse]   step %4d/%d  t=%.4g s  picard=%d  cg=%d\n",
                      step + 1, nsteps, t, std::min(picard, o.max_picard),
                      lin_iters);
        *log_ << buf;
      }
    }
  }

  if (log_ && o.verbosity >= 1) {
    for (int s = 0; s < ns; ++s) {
      double mass = 0, peak = 0;
      for (int i = 0; i < nc; ++i) {
        mass += (*fields[s].conc)[i] * mesh_.cell_vol[i];
        peak = std::max(peak, (*fields[s].conc)[i]);
      }
      char buf[200];
      if (mass0[s] > 0) {
        std::snprintf(buf, sizeof(buf),
                      "[diffuse]   %s: peak=%.4g cm^-3, dose change=%+.3g%%\n",
                      fields[s].dopant->symbol.c_str(), peak,
                      (mass - mass0[s]) / mass0[s] * 100.0);
      } else {
        std::snprintf(buf, sizeof(buf),
                      "[diffuse]   %s: peak=%.4g cm^-3, integral=%.4g atoms\n",
                      fields[s].dopant->symbol.c_str(), peak, mass);
      }
      *log_ << buf;
    }
  }
}

}  // namespace cp
