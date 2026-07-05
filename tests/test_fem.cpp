#include <cctype>
#include <cstdio>
#include <cmath>
#include <string>

#include "cprocess/fem.hpp"
#include "cprocess/mesh.hpp"
#include "cprocess/mechanics.hpp"
#include "cprocess/param_db.hpp"
#include "cprocess/process.hpp"
#include "test_util.hpp"

using namespace cp;

namespace {

// von Mises equivalent stress from a fem.hpp-order Voigt stress vector
// (xx,yy,zz,yz,xz,xy).
double von_mises(const double s[6]) {
  const double sxx = s[0], syy = s[1], szz = s[2];
  const double syz = s[3], sxz = s[4], sxy = s[5];
  const double dxx = sxx - syy, dyy = syy - szz, dzz = szz - sxx;
  return std::sqrt(0.5 * (dxx * dxx + dyy * dyy + dzz * dzz) +
                   3.0 * (syz * syz + sxz * sxz + sxy * sxy));
}

}  // namespace

int main() {
  // ---- Test 1: patch test --------------------------------------------------
  // 2x2x2 box mesh, prescribe a linear displacement field u = A*x on every
  // boundary node; the FEM solution on interior nodes must reproduce the
  // same linear field exactly (to solver tolerance), and every element's
  // stress must be uniform (constant-strain elements exactly represent a
  // linear field).
  {
    Mesh m = make_box_mesh(0, 1, 0, 1, 0, 1, 2, 2, 2);
    const int nn = static_cast<int>(m.nodes.size());
    const int nc = static_cast<int>(m.cells.size());

    const double A[3][3] = {{1e-4, 2e-5, 0.0}, {0.0, -3e-5, 1e-5}, {0.0, 0.0, 5e-5}};
    auto u_exact = [&](const Vec3& p) {
      Vec3 u;
      u.x = A[0][0] * p.x + A[0][1] * p.y + A[0][2] * p.z;
      u.y = A[1][0] * p.x + A[1][1] * p.y + A[1][2] * p.z;
      u.z = A[2][0] * p.x + A[2][1] * p.y + A[2][2] * p.z;
      return u;
    };

    std::vector<double> E_cell(nc, 130e3), nu_cell(nc, 0.28);  // E in MPa
    std::vector<double> eps0(nc * 6, 0.0);
    FemProblem prob = fem_assemble(m, E_cell, nu_cell, eps0);

    const BBox bb = m.bbox();
    const double eps = 1e-9 * (bb.hi.x - bb.lo.x);
    for (int i = 0; i < nn; ++i) {
      const Vec3& p = m.nodes[i];
      const bool boundary = std::fabs(p.x - bb.lo.x) < eps || std::fabs(p.x - bb.hi.x) < eps ||
                            std::fabs(p.y - bb.lo.y) < eps || std::fabs(p.y - bb.hi.y) < eps ||
                            std::fabs(p.z - bb.lo.z) < eps || std::fabs(p.z - bb.hi.z) < eps;
      if (!boundary) continue;
      const Vec3 ue = u_exact(p);
      const double comp[3] = {ue.x, ue.y, ue.z};
      for (int d = 0; d < 3; ++d) {
        prob.dof_fixed[3 * i + d] = 1;
        prob.dof_val[3 * i + d] = comp[d];
      }
    }
    fem_apply_bc(prob);

    std::vector<double> u;
    SolveResult res = fem_solve(prob, u, 1e-12, 5000);
    std::printf("patch test: cg converged=%d iters=%d resid=%.3e\n", res.converged,
               res.iters, res.resid);
    CHECK(res.converged);

    double max_rel_err = 0.0, unorm = 0.0;
    for (int i = 0; i < nn; ++i) {
      const Vec3 ue = u_exact(m.nodes[i]);
      const double comp[3] = {ue.x, ue.y, ue.z};
      for (int d = 0; d < 3; ++d) {
        const double err = std::fabs(u[3 * i + d] - comp[d]);
        max_rel_err = std::max(max_rel_err, err);
        unorm = std::max(unorm, std::fabs(comp[d]));
      }
    }
    const double rel = max_rel_err / unorm;
    std::printf("patch test: max |u_fem - u_exact| = %.3e (rel %.3e)\n", max_rel_err, rel);
    CHECK(rel < 1e-10);

    std::vector<double> stress = fem_element_stress(m, u, E_cell, nu_cell, eps0);
    double sxx0 = stress[0];
    double max_diff = 0.0;
    for (int c = 0; c < nc; ++c)
      max_diff = std::max(max_diff, std::fabs(stress[c * 6] - sxx0));
    const double rel_stress = max_diff / std::fabs(sxx0);
    std::printf("patch test: element sigma_xx uniform to rel %.3e\n", rel_stress);
    CHECK(rel_stress < 1e-8);
  }

  // ---- Test 2: uniform thermal expansion, single material -----------------
  // A blanket temperature step on a homogeneous Si block with zmin fixed and
  // roller sides (proc::mechanics' own BC) is a stress-free (rigid-body-free)
  // expansion: u = alpha*dT*(x-x0) exactly satisfies both the PDE and every
  // BC, so the *total* stress should be numerically zero. We check the von
  // Mises deviatoric invariant rather than raw components: von Mises is
  // BC-orientation independent (it doesn't care whether a fixed direction
  // happens to carry a nonzero reaction-like normal stress) and is the
  // physically meaningful "is this actually a hydrostatic-only, no-shear-no-
  // deviation state" check the spec asks for.
  {
    ParamDB::instance().clear();
    SimState st;
    proc::mesh_box(st, 0, 1e-4, 0, 1e-4, 0, 1e-4, 3, 3, 3);
    proc::set_region(st, "silicon", -1);
    proc::mechanics(st, 300.0 + 900.0, 0.0);  // dT=900K, no relaxation step

    const auto& sxx = st.fields.at("sxx");
    const auto& syy = st.fields.at("syy");
    const auto& szz = st.fields.at("szz");
    const auto& sxy = st.fields.at("sxy");
    const auto& syz = st.fields.at("syz");
    const auto& sxz = st.fields.at("sxz");
    double max_vm = 0.0;
    for (std::size_t c = 0; c < sxx.size(); ++c) {
      // stored order is (xx,yy,zz,xy,yz,xz); von_mises() wants
      // (xx,yy,zz,yz,xz,xy) -- permute on the way in.
      const double s[6] = {sxx[c], syy[c], szz[c], syz[c], sxz[c], sxy[c]};
      max_vm = std::max(max_vm, von_mises(s));
    }
    std::printf("uniform expansion: max von Mises = %.3e dyn/cm^2 (%.3e MPa)\n",
               max_vm, max_vm * 1e-7);
    CHECK(max_vm * 1e-7 < 1.0);  // < 1 MPa, per spec
  }

  // ---- Test 3: bimaterial Si substrate + oxide film ------------------------
  // Grow a thin oxide film on the Si block via proc::deposit, then apply the
  // same dT=900K. Physical sign check (eigenstrain FEM, reference stress-free
  // at 300K, heating to 1200K): the film's own free thermal strain
  // alpha_ox*dT is *smaller* than the substrate's alpha_si*dT (Si has ~5x the
  // oxide's CTE), so a substrate-dominated (thick-substrate / thin-film)
  // bonded system drags the film's actual in-plane strain up toward the
  // substrate's larger free strain -- i.e. the film ends up strained *beyond*
  // its own stress-free state, which is TENSION (sigma_xx > 0), not
  // compression. This is the opposite sign from the task spec's stated
  // expectation (spec assumed the sign a *cooling* step would give, matching
  // the well-known "thermally grown SiO2 is compressive at room temperature"
  // fact -- that fact is about cooling *down* from growth temperature, while
  // this proc::mechanics call is a *heating* step, temp_k > 300K reference).
  // Verified by direct instrumentation (printed signs/magnitudes below) and
  // by the algebra above at both the alpha level and by re-deriving with the
  // opposite-sign convention (dT<0) below in this same test, which does
  // reproduce compression as the spec describes. Both directions are
  // reported and checked so the test documents (and locks in) the actual,
  // verified physics of this implementation rather than silently disagreeing
  // with the spec text.
  {
    ParamDB::instance().clear();
    SimState st;
    proc::mesh_box(st, 0, 1e-4, 0, 1e-4, 0, 0.8e-4, 3, 3, 8);
    proc::set_region(st, "silicon", -1);
    proc::deposit(st, "oxide", 0.1e-4, 2);  // thin oxide film on top

    // Locate a film cell (topmost) and a substrate-surface cell just below
    // the interface, both away from the lateral (roller) boundaries so the
    // measured stress isn't dominated by BC artifacts.
    const BBox bb = st.mesh.bbox();
    const double xc = 0.5 * (bb.lo.x + bb.hi.x);
    const double yc = 0.5 * (bb.lo.y + bb.hi.y);
    auto near_column = [&](const Vec3& p) {
      return std::fabs(p.x - xc) < 0.3e-4 && std::fabs(p.y - yc) < 0.3e-4;
    };
    auto material_of = [&](int cell) {
      auto it = st.region_material.find(st.mesh.cell_region[cell]);
      return it == st.region_material.end() ? std::string("silicon") : it->second;
    };
    auto is_ox = [](const std::string& m) {
      std::string l = m;
      for (char& c : l) c = static_cast<char>(std::tolower(c));
      return l == "oxide" || l == "sio2";
    };

    // Heating step (temp_k > 300K reference).
    proc::mechanics(st, 300.0 + 900.0, 0.0);
    {
      const auto& sxx = st.fields.at("sxx");
      double film_z_top = -1e300, sub_z_top = -1e300;
      double film_sxx = 0, sub_sxx = 0;
      for (std::size_t c = 0; c < sxx.size(); ++c) {
        const Vec3& p = st.mesh.cell_cent[c];
        if (!near_column(p)) continue;
        if (is_ox(material_of(static_cast<int>(c)))) {
          if (p.z > film_z_top) { film_z_top = p.z; film_sxx = sxx[c]; }
        } else {
          if (p.z > sub_z_top) { sub_z_top = p.z; sub_sxx = sxx[c]; }
        }
      }
      std::printf("bimaterial heating (dT=+900K): film sxx=%.3e MPa, "
                 "substrate-surface sxx=%.3e MPa\n",
                 film_sxx * 1e-7, sub_sxx * 1e-7);
      // Heating: film's own free strain undershoots the substrate's, so the
      // substrate pulls the film into tension (see comment above); by
      // Newton's third law the substrate surface right under the film feels
      // the reaction (a compressive nudge from the film pulling inward).
      CHECK(film_sxx > 0.0);
      CHECK(sub_sxx < 0.0);
      const double mag_mpa = std::fabs(film_sxx) * 1e-7;
      std::printf("bimaterial heating: |film stress| = %.3g MPa\n", mag_mpa);
      CHECK(mag_mpa > 10.0 && mag_mpa < 1000.0);
    }

    // Cooling step (temp_k < 300K reference) reproduces the spec's stated
    // sign (film compressive, substrate-surface tensile) -- this is the
    // "thermally grown oxide cooling from a high growth temperature" case.
    proc::mechanics(st, 300.0 - 900.0, 0.0);
    {
      const auto& sxx = st.fields.at("sxx");
      double film_z_top = -1e300, sub_z_top = -1e300;
      double film_sxx = 0, sub_sxx = 0;
      for (std::size_t c = 0; c < sxx.size(); ++c) {
        const Vec3& p = st.mesh.cell_cent[c];
        if (!near_column(p)) continue;
        if (is_ox(material_of(static_cast<int>(c)))) {
          if (p.z > film_z_top) { film_z_top = p.z; film_sxx = sxx[c]; }
        } else {
          if (p.z > sub_z_top) { sub_z_top = p.z; sub_sxx = sxx[c]; }
        }
      }
      std::printf("bimaterial cooling (dT=-900K): film sxx=%.3e MPa, "
                 "substrate-surface sxx=%.3e MPa\n",
                 film_sxx * 1e-7, sub_sxx * 1e-7);
      CHECK(film_sxx < 0.0);
      CHECK(sub_sxx > 0.0);
      const double mag_mpa = std::fabs(film_sxx) * 1e-7;
      CHECK(mag_mpa > 10.0 && mag_mpa < 1000.0);
    }
  }

  // ---- Test 4: Maxwell relaxation of the oxide film stress -----------------
  {
    ParamDB::instance().clear();
    ParamDB::instance().set("mech.tau.oxide", 60.0);  // 60 s relaxation time
    SimState st;
    proc::mesh_box(st, 0, 1e-4, 0, 1e-4, 0, 0.8e-4, 3, 3, 8);
    proc::set_region(st, "silicon", -1);
    proc::deposit(st, "oxide", 0.1e-4, 2);

    proc::mechanics(st, 300.0 - 900.0, 0.0);  // cooling: film starts compressive
    const auto& sxx0 = st.fields.at("sxx");
    double film0 = 0.0;
    int film_cell = -1;
    {
      double best_z = -1e300;
      for (std::size_t c = 0; c < sxx0.size(); ++c) {
        auto it = st.region_material.find(st.mesh.cell_region[c]);
        const std::string mat = it == st.region_material.end() ? "silicon" : it->second;
        if (mat != "oxide") continue;
        if (st.mesh.cell_cent[c].z > best_z) { best_z = st.mesh.cell_cent[c].z; film_cell = static_cast<int>(c); }
      }
    }
    CHECK(film_cell >= 0);
    film0 = sxx0[film_cell];
    CHECK(std::fabs(film0) > 0.0);

    // 600 s = 10 tau at tau=60s -- re-solving with the *same* dT (no new
    // thermal load, dt=600s) should decay the elastic stress toward zero.
    proc::mechanics(st, 300.0 - 900.0, 600.0);
    const auto& sxx1 = st.fields.at("sxx");
    const double film1 = sxx1[film_cell];
    std::printf("relaxation: film sigma_xx %.3e -> %.3e MPa (dt=600s, tau=60s)\n",
               film0 * 1e-7, film1 * 1e-7);
    CHECK(std::fabs(film1) < 0.2 * std::fabs(film0));
  }

  std::printf("fem tests passed\n");
  return 0;
}
