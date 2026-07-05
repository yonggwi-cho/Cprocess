#include "cprocess/oxidant_solver.hpp"

#include <algorithm>
#include <cmath>
#include <ostream>

#include "cprocess/sparse.hpp"

namespace cp {

// Implementation note (P2-4 Stage A): a plain per-tet two-point flux badly
// overestimates the oxide's series resistance on a Kuhn-tet box mesh --
// every hex is split into 6 tets sharing its main space diagonal, and the
// internal faces between those 6 tets are severely non-orthogonal (S is not
// remotely aligned with the owner->neigh line). A deferred (explicit)
// non-orthogonal correction, Picard-iterated the way DiffusionSolver does
// for the transient dopant solve, was tried first but is unconditionally
// unstable here: the transient solve's V/dt diagonal damps that iteration,
// and a pure steady Laplace system has no such term, so the correction
// amplifies along any multi-hex chain (verified empirically while
// developing this solver: stable for 1-2 layers, diverges by 5-10).
//
// Fix: since every mesh this solver is used on (make_box_mesh output, per
// the P2-4 spec's own mesh model) has a fixed, regular structure -- cells
// are emitted in groups of 6 per hex, in hex-scan order -- solve at
// hex-aggregate resolution instead of raw-tet resolution. Faces *within* a
// hex (the problematic diagonal ones) are never part of the hex-level
// system at all; faces *between* hexes are the true axis-aligned quad
// faces (split into 2 triangles by the same diagonal, but each triangle's
// own centroid-to-centroid line is still close enough to its area vector
// that summing the 2 triangles' transmissibilities in parallel reproduces
// the exact axis-aligned quad conductance, no clamping ever triggered).
// This is exact for the intended regular-grid use case and is documented
// here as a deliberate scope decision rather than a fully general
// unstructured-mesh solver.
namespace {

struct IFace { int oh, nh; double t; };            // inter-hex internal face
struct DFace { int oh; double tb; };                // Dirichlet gas (zmax)
struct RFace { int oh; double t_robin, area; int face; };  // Robin, per tet-face

int hex_of(int cell) { return cell / 6; }

}  // namespace

OxidantResult solve_oxidant(const Mesh& m, const std::vector<char>& oxide_mask,
                            const std::vector<char>& si_mask, double d_ox,
                            double ks, double c_gas, std::ostream* log) {
  std::vector<double> d_cell(m.cells.size(), d_ox);
  return solve_oxidant(m, d_cell, oxide_mask, si_mask, ks, c_gas, log);
}

OxidantResult solve_oxidant(const Mesh& m, const std::vector<double>& d_cell,
                            const std::vector<char>& active_mask,
                            const std::vector<char>& si_mask, double ks,
                            double c_gas, std::ostream* log) {
  const int nc = static_cast<int>(m.cells.size());
  const int nf = static_cast<int>(m.faces.size());
  OxidantResult res;
  res.conc.assign(nc, 0.0);
  if (nc == 0) return res;
  if (nc % 6 != 0)
    throw std::runtime_error("solve_oxidant: expects a make_box_mesh-style "
                              "Kuhn-tet mesh (cell count not a multiple of 6)");
  const int nh = nc / 6;

  // Per-hex aggregates: active/D come from any active tet in the group (by
  // construction all 6 share the same region/material within a hex).
  std::vector<char> active_hex(nh, 0);
  std::vector<double> d_hex(nh, 0.0);
  for (int i = 0; i < nc; ++i) {
    const int h = hex_of(i);
    if (active_mask[i]) { active_hex[h] = 1; d_hex[h] = d_cell[i]; }
  }
  // Hex centroid: average of its 6 tet centroids (exact center for the
  // symmetric Kuhn split of a parallelepiped hex).
  std::vector<Vec3> hex_cent(nh, Vec3{});
  for (int i = 0; i < nc; ++i) hex_cent[hex_of(i)] += (1.0 / 6.0) * m.cell_cent[i];

  std::vector<IFace> ifaces;
  std::vector<DFace> dfaces;
  std::vector<RFace> rfaces;

  for (int fi = 0; fi < nf; ++fi) {
    const Face& f = m.faces[fi];
    const int o = f.owner, n = f.neigh;
    const double area = norm(f.S);
    if (area <= 0) continue;
    const int ho = hex_of(o);
    const bool ao = active_mask[o];
    const bool an = (n >= 0) && active_mask[n];

    if (n >= 0 && hex_of(n) == ho) continue;  // intra-hex diagonal face: skip

    if (ao && an) {
      const int hn = hex_of(n);
      const double d0 = d_cell[o], d1 = d_cell[n];
      if (d0 <= 0 || d1 <= 0) continue;
      const Vec3 dvec = hex_cent[hn] - hex_cent[ho];
      const double sn = area, dn = norm(dvec);
      const double sd = std::max(dot(f.S, dvec), 1e-12 * sn * dn + 1e-300);
      const double g = dot(f.S, f.S) / sd;
      const double dh = (d0 > 0 && d1 > 0) ? (2.0 * d0 * d1 / (d0 + d1)) : 0.0;
      ifaces.push_back({ho, hn, dh * g});
    } else if (ao && n < 0) {
      const std::string pname =
          (f.patch >= 0 && f.patch < static_cast<int>(m.patch_names.size()))
              ? m.patch_names[f.patch] : "";
      if (pname == "zmax") {
        const Vec3 d = f.c - hex_cent[ho];
        const double sn = area, dn = norm(d);
        const double sd = std::max(dot(f.S, d), 1e-12 * sn * dn + 1e-300);
        const double gb = dot(f.S, f.S) / sd;
        dfaces.push_back({ho, d_cell[o] * gb});
      }
    } else if (ao && n >= 0 && !an && si_mask[n]) {
      const Vec3 d = f.c - hex_cent[ho];
      const double sn = area, dn = norm(d);
      const double sd = std::max(dot(f.S, d), 1e-12 * sn * dn + 1e-300);
      const double gb = dot(f.S, f.S) / sd;
      const double t_diff = d_cell[o] * gb;
      if (t_diff > 0) {
        const double t_robin = 1.0 / (1.0 / t_diff + 1.0 / (ks * area));
        rfaces.push_back({ho, t_robin, area, fi});
      }
    } else if (an && n >= 0 && si_mask[o] && !ao) {
      const int hn = hex_of(n);
      const Vec3 d = f.c - hex_cent[hn];
      const double sn = area, dn = norm(d);
      const double sd = std::max(-dot(f.S, d), 1e-12 * sn * dn + 1e-300);
      const double gb = dot(f.S, f.S) / sd;
      const double t_diff = d_cell[n] * gb;
      if (t_diff > 0) {
        const double t_robin = 1.0 / (1.0 / t_diff + 1.0 / (ks * area));
        rfaces.push_back({hn, t_robin, area, fi});
      }
    }
  }

  // Hex-level CSR.
  CSR A;
  A.n = nh;
  std::vector<std::vector<int>> nb(nh);
  for (int h = 0; h < nh; ++h) nb[h].push_back(h);
  for (const auto& fc : ifaces) { nb[fc.oh].push_back(fc.nh); nb[fc.nh].push_back(fc.oh); }
  A.ptr.assign(nh + 1, 0);
  for (int h = 0; h < nh; ++h) {
    auto& v = nb[h];
    std::sort(v.begin(), v.end());
    v.erase(std::unique(v.begin(), v.end()), v.end());
    A.ptr[h + 1] = A.ptr[h] + static_cast<int>(v.size());
  }
  A.col.resize(A.ptr[nh]);
  for (int h = 0; h < nh; ++h) std::copy(nb[h].begin(), nb[h].end(), A.col.begin() + A.ptr[h]);
  A.val.assign(A.ptr[nh], 0.0);
  std::vector<int> diag(nh);
  for (int h = 0; h < nh; ++h) diag[h] = A.find(h, h);

  std::vector<double> rhs(nh, 0.0);
  for (int h = 0; h < nh; ++h) if (!active_hex[h]) A.val[diag[h]] = 1.0;
  for (const auto& fc : ifaces) {
    A.val[diag[fc.oh]] += fc.t;
    A.val[diag[fc.nh]] += fc.t;
    A.val[A.find(fc.oh, fc.nh)] -= fc.t;
    A.val[A.find(fc.nh, fc.oh)] -= fc.t;
  }
  for (const auto& fc : dfaces) { A.val[diag[fc.oh]] += fc.tb; rhs[fc.oh] += fc.tb * c_gas; }
  for (const auto& fc : rfaces) A.val[diag[fc.oh]] += fc.t_robin;

  std::vector<double> x(nh, 0.0);
  const SolveResult sr = cg_ilu0(A, rhs, x, 1e-10, 2000);
  if (log && !sr.converged)
    *log << "[oxidant] warning: cg_ilu0 did not converge to rtol (resid="
         << sr.resid << ", iters=" << sr.iters << ")\n";

  for (int i = 0; i < nc; ++i) res.conc[i] = active_mask[i] ? x[hex_of(i)] : 0.0;

  constexpr double kNOx = 2.25e22;  // SiO2 molecular density [cm^-3]
  res.iface_vn.reserve(rfaces.size());
  for (const auto& fc : rfaces) {
    const double F = fc.t_robin * x[fc.oh];
    const double vn = F / (fc.area * kNOx);
    res.iface_vn.emplace_back(fc.face, vn);
  }
  return res;
}

}  // namespace cp
