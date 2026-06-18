// Physical photoresist masking for MC implantation.
//
// A box is split into a top resist layer and a Si substrate below it. Half of
// the lateral area (x < 0.5) is covered by resist (material 1); the open half
// (x >= 0.5) is silicon all the way up (material 0). A full-surface implant
// must then deposit far more dopant into the open silicon than under the
// resist, and the resist must stop a substantial fraction of ions.

#include <cstdio>
#include <vector>

#include "cprocess/mc_implant.hpp"
#include "cprocess/mesh.hpp"
#include "test_util.hpp"

using namespace cp;

int main() {
  // 1 um x 1 um x 0.6 um; top 0.2 um is the resist layer over x < 0.5 um.
  const double Lx = 1e-4, Ly = 1e-4, Lz = 0.6e-4;
  const double z_resist_bot = 0.2e-4;  // resist occupies z in [0.2, 0.6] um (0.4 um thick)
  const double x_split = 0.5e-4;
  Mesh m = make_box_mesh(0, Lx, 0, Ly, 0, Lz, 10, 6, 12);
  const int nc = static_cast<int>(m.cells.size());

  // Material map: 0 = Si (substrate + open column), 1 = resist.
  std::vector<int> cell_material(nc, 0);
  std::vector<char> si_mask(nc, 1);
  for (int ci = 0; ci < nc; ++ci) {
    const Vec3& c = m.cell_cent[ci];
    if (c.z > z_resist_bot && c.x < x_split) {
      cell_material[ci] = 1;  // resist
      si_mask[ci] = 0;        // dopant resting in resist is discarded
    }
  }

  const Dopant* B = find_dopant("B");
  CHECK(B != nullptr);

  McImplantParams p;
  p.dopant         = B;
  p.dose           = 1e15;
  p.energy_kev     = 30;     // low energy: resist (0.2 um) blocks most ions
  p.ions           = 60000;
  p.channeling     = false;  // deterministic; isolate the masking effect
  p.seed           = 99;
  p.material_table = {target_silicon(), target_photoresist()};
  p.cell_material  = &cell_material;

  std::vector<double> conc(nc, 0.0);
  const auto st = apply_mc_implant(m, si_mask, p, conc, nullptr);
  std::printf("resist run: deposited=%lld in_resist=%lld backscattered=%lld\n",
              st.deposited, st.in_mask, st.backscattered);

  // Integrate dopant that reached silicon (si_mask cells) in each lateral half.
  double dose_open = 0, dose_masked = 0;
  for (int ci = 0; ci < nc; ++ci) {
    if (!si_mask[ci]) continue;  // silicon only
    const double n = conc[ci] * m.cell_vol[ci];
    if (m.cell_cent[ci].x < x_split) dose_masked += n;
    else                             dose_open   += n;
  }
  std::printf("silicon dopant: open=%.3e  masked=%.3e  ratio=%.1fx\n",
              dose_open, dose_masked, dose_open / (dose_masked + 1e-30));

  // The resist must absorb a meaningful fraction of ions.
  CHECK(st.in_mask > p.ions / 20);
  // The open silicon must receive far more dopant than under the mask.
  CHECK(dose_open > 5.0 * dose_masked);

  // Control: same implant with no resist deposits everywhere roughly equally.
  {
    std::vector<int> all_si(nc, 0);
    std::vector<char> all_mask(nc, 1);
    McImplantParams p2 = p;
    p2.material_table = {target_silicon()};
    p2.cell_material  = &all_si;
    std::vector<double> conc2(nc, 0.0);
    apply_mc_implant(m, all_mask, p2, conc2, nullptr);
    double left = 0, right = 0;
    for (int ci = 0; ci < nc; ++ci) {
      const double n = conc2[ci] * m.cell_vol[ci];
      if (m.cell_cent[ci].x < x_split) left += n; else right += n;
    }
    std::printf("control (no resist): left=%.3e right=%.3e ratio=%.2f\n",
                left, right, left / (right + 1e-30));
    // Without a mask the two halves are within ~30%.
    CHECK(left > 0.7 * right && left < 1.3 * right);
  }

  std::printf("resist tests passed\n");
  return 0;
}
