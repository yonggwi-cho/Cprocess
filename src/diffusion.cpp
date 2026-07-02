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
                                 std::ostream* log, std::vector<int> cell_mat)
    : mesh_(mesh), mask_(std::move(solve_mask)), log_(log) {
  const int nc = static_cast<int>(mesh_.cells.size());
  if (!cell_mat.empty()) {
    if (static_cast<int>(cell_mat.size()) != nc)
      throw std::runtime_error("diffusion: cell_mat size mismatch");
    mat_ = std::move(cell_mat);
    mask_.assign(nc, 0);
    for (int i = 0; i < nc; ++i) mask_[i] = (mat_[i] == 0 || mat_[i] == 1) ? 1 : 0;
  } else {
    if (mask_.empty()) mask_.assign(nc, 1);
    if (static_cast<int>(mask_.size()) != nc)
      throw std::runtime_error("diffusion: mask size mismatch");
    mat_.assign(nc, 0);
    for (int i = 0; i < nc; ++i) mat_[i] = mask_[i] ? 0 : 2;
  }
  build();
}

void DiffusionSolver::build() {
  const int nc = static_cast<int>(mesh_.cells.size());
  const int nf = static_cast<int>(mesh_.faces.size());
  fg_.assign(nf, FGeom{});
  clamped_faces_ = 0;
  has_segregation_ = false;

  for (int fi = 0; fi < nf; ++fi) {
    const Face& f = mesh_.faces[fi];
    const bool mo = mask_[f.owner];
    const bool mn = f.neigh >= 0 && mask_[f.neigh];
    FGeom& g = fg_[fi];
    const bool sameMat = mo && mn && mat_[f.owner] == mat_[f.neigh];
    const bool isSeg = mo && mn && mat_[f.owner] != mat_[f.neigh];
    if (sameMat) {
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
    } else if (isSeg) {
      g.kind = kSegregation;
      g.area = norm(f.S);
      has_segregation_ = true;
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
    if (fg_[fi].kind != kInternal && fg_[fi].kind != kSegregation) continue;
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
    if (fg_[fi].kind != kInternal && fg_[fi].kind != kSegregation) continue;
    const Face& f = mesh_.faces[fi];
    fslot_[fi] = {A_.find(f.owner, f.neigh), A_.find(f.neigh, f.owner)};
  }

  // Greedy coloring of active faces. Two faces conflict if they share a cell
  // (owner or neigh); giving each face the lowest color used by none of its
  // already-colored conflicting faces guarantees that same-color faces write to
  // disjoint cell rows, so parallel scatter-add needs no atomics.
  //
  // A face contributes to a cell row iff it writes that cell's diagonal/rhs:
  //   kInternal, kSegregation -> owner and neigh
  //   kBoundOwner             -> owner only
  // (kBoundNeigh / kInactive touch no active row and are skipped.)
  auto writes_cells = [&](int fi, int& c0, int& c1) {
    const Face& f = mesh_.faces[fi];
    c0 = f.owner;
    c1 = (fg_[fi].kind == kInternal || fg_[fi].kind == kSegregation) ? f.neigh : -1;
  };
  std::vector<std::vector<int>> cell_faces(nc);
  std::vector<int> active_faces;
  for (int fi = 0; fi < nf; ++fi) {
    if (fg_[fi].kind != kInternal && fg_[fi].kind != kBoundOwner &&
        fg_[fi].kind != kSegregation) continue;
    int c0, c1;
    writes_cells(fi, c0, c1);
    cell_faces[c0].push_back(fi);
    if (c1 >= 0) cell_faces[c1].push_back(fi);
    active_faces.push_back(fi);
  }
  std::vector<int> face_color(nf, -1);
  std::vector<int> forbidden;  // forbidden[color] == fi means color taken
  int ncolors = 0;
  for (int fi : active_faces) {
    int c0, c1;
    writes_cells(fi, c0, c1);
    // Mark colors already used by faces sharing either incident cell.
    for (int cc : {c0, c1}) {
      if (cc < 0) continue;
      for (int gf : cell_faces[cc]) {
        const int gc = face_color[gf];
        if (gc < 0) continue;
        if (gc >= static_cast<int>(forbidden.size())) forbidden.resize(gc + 1, -1);
        forbidden[gc] = fi;
      }
    }
    int color = 0;
    while (color < static_cast<int>(forbidden.size()) && forbidden[color] == fi)
      ++color;
    face_color[fi] = color;
    ncolors = std::max(ncolors, color + 1);
  }
  face_colors_.assign(ncolors, {});
  for (int fi : active_faces) face_colors_[face_color[fi]].push_back(fi);

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

  // Face-colored scatter so owner/neigh accumulations never race in parallel.
  // kBoundNeigh faces write only the neigh cell and are excluded from the
  // coloring (they touch no "active row" for assembly), so handle them in a
  // separate serial pass — they are few (material interfaces) and cheap.
  for (const auto& color : face_colors_) {
    const int mfc = static_cast<int>(color.size());
#pragma omp parallel for schedule(static)
    for (int t = 0; t < mfc; ++t) {
      const int fi = color[t];
      const Face& f = mesh_.faces[fi];
      if (fg_[fi].kind == kInternal) {
        const Vec3 d = mesh_.cell_cent[f.neigh] - mesh_.cell_cent[f.owner];
        const double dc = c[f.neigh] - c[f.owner];
        add_row(f.owner, d, dc);
        add_row(f.neigh, -1.0 * d, -dc);
      } else if (fg_[fi].kind == kBoundOwner) {
        const Vec3 d = f.c - mesh_.cell_cent[f.owner];
        if (!std::isnan(bcface[fi])) {
          add_row(f.owner, d, bcface[fi] - c[f.owner]);
        } else {
          const Vec3 nh = (1.0 / norm(f.S)) * f.S;
          add_row(f.owner, dot(d, nh) * nh, 0.0);
        }
      }
    }
  }
  for (std::size_t fi = 0; fi < mesh_.faces.size(); ++fi) {
    if (fg_[fi].kind != kBoundNeigh) continue;
    const Face& f = mesh_.faces[fi];
    const Vec3 d = f.c - mesh_.cell_cent[f.neigh];
    const Vec3 nh = (1.0 / norm(f.S)) * f.S;
    add_row(f.neigh, dot(d, nh) * nh, 0.0);
  }

#pragma omp parallel for schedule(static)
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

void DiffusionSolver::assemble(const std::vector<double>& dcell,
                               const std::vector<double>& cold,
                               const std::vector<double>& bcface,
                               const std::vector<double>& cgrad, double dt,
                               double reaction, bool nonortho,
                               double h_seg, double m_seg,
                               std::vector<double>& rhs, std::vector<Vec3>& grad) {
  const int nc = static_cast<int>(mesh_.cells.size());
  gradients(cgrad, bcface, grad);

  std::fill(A_.val.begin(), A_.val.end(), 0.0);
#pragma omp parallel for schedule(static)
  for (int i = 0; i < nc; ++i) {
    if (!mask_[i]) {
      A_.val[diag_[i]] = 1.0;
      rhs[i] = cgrad[i];  // frozen cell holds its value
    } else {
      const double vd = mesh_.cell_vol[i] / dt;
      A_.val[diag_[i]] = vd + reaction * mesh_.cell_vol[i];
      rhs[i] = vd * cold[i];
    }
  }
  // Scatter face fluxes color by color: within a color no two faces share a
  // cell, so the diagonal/rhs updates are race-free. The off-diagonal slots
  // fslot_[fi][0/1] are unique to a face and never conflict.
  for (const auto& color : face_colors_) {
    const int mfc = static_cast<int>(color.size());
#pragma omp parallel for schedule(static)
    for (int t = 0; t < mfc; ++t) {
      const int fi = color[t];
      const FGeom& g = fg_[fi];
      const Face& f = mesh_.faces[fi];
      if (g.kind == kInternal) {
        const double dP = dcell[f.owner], dN = dcell[f.neigh];
        if (dP <= 0 || dN <= 0) continue;
        const double dh = (g.delP + g.delN) / (g.delP / dP + g.delN / dN);
        const double tf = dh * g.g;
        A_.val[diag_[f.owner]] += tf;
        A_.val[diag_[f.neigh]] += tf;
        A_.val[fslot_[fi][0]] -= tf;
        A_.val[fslot_[fi][1]] -= tf;
        if (nonortho) {
          const Vec3 gf = g.wP * grad[f.owner] + (1.0 - g.wP) * grad[f.neigh];
          const double corr = dh * dot(gf, g.k);
          rhs[f.owner] += corr;
          rhs[f.neigh] -= corr;
        }
      } else if (g.kind == kBoundOwner && !std::isnan(bcface[fi])) {
        const double tb = dcell[f.owner] * g.gb;
        A_.val[diag_[f.owner]] += tb;
        rhs[f.owner] += tb * bcface[fi];
      } else if (g.kind == kSegregation && h_seg > 0) {
        // owner/neigh may be Si or oxide in either order; identify iS (Si)
        // and iO (oxide) via mat_, and pick the fslot_ slots accordingly.
        const bool ownerIsSi = mat_[f.owner] == 0;
        const int iS = ownerIsSi ? f.owner : f.neigh;
        const int iO = ownerIsSi ? f.neigh : f.owner;
        const double hA = h_seg * g.area;
        // fslot_[fi][0] is the (owner,neigh) slot, [1] is (neigh,owner).
        const int s_S_to_O = ownerIsSi ? fslot_[fi][0] : fslot_[fi][1];
        const int s_O_to_S = ownerIsSi ? fslot_[fi][1] : fslot_[fi][0];
        A_.val[diag_[iS]]  += hA;
        A_.val[s_S_to_O]   -= hA * m_seg;
        A_.val[diag_[iO]]  += hA * m_seg;
        A_.val[s_O_to_S]   -= hA;
      }
    }
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
      // Oxide cells carry no free carriers relevant to this model: nni=1.
      for (int i = 0; i < nc; ++i) {
        if (mat_[i] != 0) { nni[i] = 1.0; continue; }
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
          if (mat_[i] == 1) { dcell[s][i] = oxide_diffusivity(dp, o.temp); continue; }
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
        const Dopant& dp = *fields[s].dopant;
        std::vector<double>& c = *fields[s].conc;
        const double h_seg = has_segregation_ ? segregation_h(dp, o.temp) : 0.0;
        const double m_seg = segregation_m(dp, o.temp);
        assemble(dcell[s], cold[s], bcface[s], c, dt, 0.0, o.nonortho, h_seg,
                 m_seg, rhs, grad);

        x = c;  // warm start
        // ILU(0)-preconditioned CG (parallel, level-scheduled triangular solves)
        // with a BiCGSTAB fallback for the occasional non-SPD assembly. The
        // segregation cross-terms make the matrix non-symmetric when m != 1,
        // so skip straight to BiCGSTAB in that case.
        SolveResult sr;
        if (has_segregation_ && h_seg > 0) {
          sr = bicgstab_ilu0(A_, rhs, x, o.lin_rtol, o.lin_maxit);
        } else {
          sr = cg_ilu0(A_, rhs, x, o.lin_rtol, o.lin_maxit);
        }
        if (!sr.converged) {
          x = c;
          sr = bicgstab_ilu0(A_, rhs, x, o.lin_rtol, o.lin_maxit);
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

void DiffusionSolver::run_ted(std::vector<SpeciesField>& fields,
                              std::vector<double>& psi,
                              const std::vector<DirichletBC>& bcs,
                              const DiffuseOpts& o) {
  const int nc = static_cast<int>(mesh_.cells.size());
  const int nf = static_cast<int>(mesh_.faces.size());
  const int ns = static_cast<int>(fields.size());
  for (auto& sf : fields) {
    if (!sf.dopant || !sf.conc) throw std::runtime_error("ted: bad field");
    sf.conc->resize(nc, 0.0);
  }
  psi.resize(nc, 0.0);
  std::sort(fields.begin(), fields.end(),
            [](const SpeciesField& a, const SpeciesField& b) {
              return a.dopant->symbol < b.dopant->symbol;
            });

  const double dt0 = (o.dt > 0) ? std::min(o.dt, o.time) : o.time / 50.0;
  const int nsteps = static_cast<int>(std::ceil(o.time / dt0 - 1e-12));
  const double ni = ni_si(o.temp);

  // Interstitial parameters (isothermal, spatially uniform).
  const double cstar = interstitial_cstar(o.temp);
  const double d_I = interstitial_diffusivity(o.temp);
  const double k_rec = interstitial_recomb_rate(o.temp);

  // Per-species Dirichlet value per boundary face (NaN = none).
  std::vector<std::vector<double>> bcface(ns, std::vector<double>(nf, kNaN));
  for (int s = 0; s < ns; ++s)
    for (const auto& bc : bcs) {
      if (bc.species != fields[s].dopant->symbol || bc.patch < 0) continue;
      for (int fi = 0; fi < nf; ++fi)
        if (mesh_.faces[fi].patch == bc.patch && fg_[fi].kind == kBoundOwner &&
            mesh_.faces[fi].neigh < 0)
          bcface[s][fi] = bc.conc;
    }

  // Interstitial surface sink: excess psi = 0 at the top (zmax) surface, where
  // interstitials recombine at the surface. This drives the transient decay and
  // TED's characteristic depth asymmetry.
  const int ztop = mesh_.find_patch("zmax");
  std::vector<double> psi_bcface(nf, kNaN);
  for (int fi = 0; fi < nf; ++fi)
    if (mesh_.faces[fi].patch == ztop && fg_[fi].kind == kBoundOwner &&
        mesh_.faces[fi].neigh < 0)
      psi_bcface[fi] = 0.0;

  // Interstitial diffusivity is constant in silicon, zero in frozen regions.
  std::vector<double> dI(nc, 0.0);
  for (int i = 0; i < nc; ++i) dI[i] = mask_[i] ? d_I : 0.0;

  std::vector<double> mass0(ns, 0.0);
  for (int s = 0; s < ns; ++s)
    for (int i = 0; i < nc; ++i)
      mass0[s] += (*fields[s].conc)[i] * mesh_.cell_vol[i];

  std::vector<std::vector<double>> cold(ns), dcell(ns, std::vector<double>(nc, 0.0));
  std::vector<double> psi_old(nc), nni(nc, 1.0), rhs(nc), x;
  std::vector<Vec3> grad;

  double t = 0, s_peak0 = 0;
  for (int i = 0; i < nc; ++i)
    s_peak0 = std::max(s_peak0, 1.0 + psi[i] / cstar);

  for (int step = 0; step < nsteps; ++step) {
    const double dt = std::min(dt0, o.time - t);

    // ── 1. Advance the interstitial excess one implicit step (linear). ──
    psi_old = psi;
    assemble(dI, psi_old, psi_bcface, psi_old, dt, k_rec, o.nonortho, 0.0, 0.0,
             rhs, grad);
    x = psi;
    SolveResult sp = cg_ilu0(A_, rhs, x, o.lin_rtol, o.lin_maxit);
    if (!sp.converged) { x = psi; bicgstab_ilu0(A_, rhs, x, o.lin_rtol, o.lin_maxit); }
    psi = x;
    for (double& v : psi) v = std::max(v, 0.0);

    // Supersaturation S = 1 + psi/C_I* (>= 1) from the updated interstitials.
    // The "+1" model over-counts free interstitials: most cluster into {311}
    // defects that buffer the free concentration. Lacking an explicit cluster
    // model, cap S at kSmax to keep the free supersaturation physical (~1e3).
    constexpr double kSmax = 3.0e3;
    double smax = 0;
    std::vector<double> S(nc, 1.0);
    for (int i = 0; i < nc; ++i) {
      S[i] = 1.0 + (mask_[i] ? std::min(psi[i] / cstar, kSmax) : 0.0);
      smax = std::max(smax, S[i]);
    }

    // ── 2. Advance dopants with interstitial-enhanced diffusivity. ──
    if (ns > 0) {
      for (int s = 0; s < ns; ++s) cold[s] = *fields[s].conc;
      double cmax = 0;
      for (int s = 0; s < ns; ++s)
        for (int i = 0; i < nc; ++i) cmax = std::max(cmax, (*fields[s].conc)[i]);
      const double cfloor = 1e-3 * std::max(cmax, 1.0);

      for (int picard = 1; picard <= o.max_picard; ++picard) {
        for (int i = 0; i < nc; ++i) {
          if (mat_[i] != 0) { nni[i] = 1.0; continue; }
          double nnet = 0;
          for (int s = 0; s < ns; ++s) {
            const double c = (*fields[s].conc)[i];
            nnet += (fields[s].dopant->type == DopType::donor) ? c : -c;
          }
          const double cc = nnet / (2.0 * ni);
          nni[i] = cc + std::sqrt(cc * cc + 1.0);
        }
        for (int s = 0; s < ns; ++s) {
          const Dopant& dp = *fields[s].dopant;
          for (int i = 0; i < nc; ++i) {
            if (!mask_[i]) { dcell[s][i] = 0; continue; }
            if (mat_[i] == 1) { dcell[s][i] = oxide_diffusivity(dp, o.temp); continue; }
            double dv = dopant_diffusivity(dp, o.temp, nni[i]);
            if (o.field_enh) {
              const bool ntype = nni[i] >= 1.0;
              if ((dp.type == DopType::donor && ntype) ||
                  (dp.type == DopType::acceptor && !ntype)) {
                const double cc = 0.5 * (nni[i] - 1.0 / nni[i]);
                dv *= 1.0 + std::fabs(cc) / std::sqrt(cc * cc + 1.0);
              }
            }
            // Pair-diffusion enhancement: (1 - fi) + fi * S.
            dv *= (1.0 - dp.fi) + dp.fi * S[i];
            dcell[s][i] = dv;
          }
        }

        double maxrel = 0;
        for (int s = 0; s < ns; ++s) {
          const Dopant& dp = *fields[s].dopant;
          std::vector<double>& c = *fields[s].conc;
          const double h_seg = has_segregation_ ? segregation_h(dp, o.temp) : 0.0;
          const double m_seg = segregation_m(dp, o.temp);
          assemble(dcell[s], cold[s], bcface[s], c, dt, 0.0, o.nonortho, h_seg,
                   m_seg, rhs, grad);
          x = c;
          SolveResult sr;
          if (has_segregation_ && h_seg > 0) {
            sr = bicgstab_ilu0(A_, rhs, x, o.lin_rtol, o.lin_maxit);
          } else {
            sr = cg_ilu0(A_, rhs, x, o.lin_rtol, o.lin_maxit);
          }
          if (!sr.converged) {
            x = c;
            sr = bicgstab_ilu0(A_, rhs, x, o.lin_rtol, o.lin_maxit);
            if (!sr.converged)
              throw std::runtime_error("ted: linear solver failed");
          }
          for (int i = 0; i < nc; ++i)
            maxrel = std::max(maxrel,
                              std::fabs(x[i] - c[i]) / (std::fabs(x[i]) + cfloor));
          c = x;
        }
        if (maxrel < o.picard_tol) break;
      }
      for (int s = 0; s < ns; ++s)
        for (double& v : *fields[s].conc) v = std::max(v, 0.0);
    }

    t += dt;
    if (o.verbosity >= 1 && log_) {
      const int every = std::max(1, nsteps / 10);
      if (step % every == 0 || step == nsteps - 1) {
        char buf[160];
        std::snprintf(buf, sizeof(buf),
                      "[ted]   step %4d/%d  t=%.4g s  Smax=%.3g\n",
                      step + 1, nsteps, t, smax);
        *log_ << buf;
      }
    }
  }

  if (log_ && o.verbosity >= 1) {
    char b0[96];
    std::snprintf(b0, sizeof(b0),
                  "[ted]   initial Smax=%.3g, C_I*=%.3g cm^-3\n", s_peak0, cstar);
    *log_ << b0;
    for (int s = 0; s < ns; ++s) {
      double mass = 0, peak = 0;
      for (int i = 0; i < nc; ++i) {
        mass += (*fields[s].conc)[i] * mesh_.cell_vol[i];
        peak = std::max(peak, (*fields[s].conc)[i]);
      }
      char buf[200];
      std::snprintf(buf, sizeof(buf),
                    "[ted]   %s: peak=%.4g cm^-3, dose change=%+.3g%%\n",
                    fields[s].dopant->symbol.c_str(), peak,
                    mass0[s] > 0 ? (mass - mass0[s]) / mass0[s] * 100.0 : 0.0);
      *log_ << buf;
    }
  }
}

}  // namespace cp
