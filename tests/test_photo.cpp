// Test for photo / mask / strip deck commands.
//
// Verifies that a physical photoresist deposited via 'photo', patterned with
// 'mask', and then used in an MC implant blocks ions under the unexposed
// region and lets them through in the opening.

#include <cstdio>
#include <sstream>

#include "cprocess/deck.hpp"
#include "cprocess/mesh.hpp"
#include "test_util.hpp"

using namespace cp;

int main() {
  // ── Deck with photo/mask/strip ─────────────────────────────────────────
  // 1 um × 1 um × 0.5 um Si substrate; 0.4 um resist; gate [0.35, 0.65 um].
  const char* deck = R"(
mesh box xmax=1um ymax=1um zmax=0.5um nx=8 ny=8 nz=4
region all material=silicon
init species=B conc=1e15
photo resist=0.4um nz=4
mask x1=0um x2=0.35um
mask x1=0.65um x2=1um
implant species=P dose=5e14 energy=30keV method=mc ions=40000 seed=17
strip
save file=/tmp/test_photo_resist.vtu
print
)";

  SimState st;
  std::istringstream in(deck);
  std::ostringstream log;
  bool ok = true;
  try {
    run_deck(in, st, log);
  } catch (const std::exception& e) {
    std::fprintf(stderr, "DECK ERROR: %s\n", e.what());
    ok = false;
  }
  std::printf("%s", log.str().c_str());
  CHECK(ok);
  CHECK(st.has_mesh);
  CHECK(!st.has_stack);  // strip was called

  // Phosphorus must be present.
  CHECK(st.fields.count("P") > 0);
  const auto& p_conc = st.fields.at("P");

  // Compute dose in open S/D region (x < 0.35 um or x > 0.65 um) vs
  // under gate (0.35 < x < 0.65).
  double dose_open = 0, dose_gate = 0;
  const int nc = static_cast<int>(st.mesh.cells.size());
  for (int ci = 0; ci < nc; ++ci) {
    const double x = st.mesh.cell_cent[ci].x;
    const double q = p_conc[ci] * st.mesh.cell_vol[ci];
    if (x < 0.35e-4 || x > 0.65e-4) dose_open += q;
    else                             dose_gate += q;
  }
  std::printf("P dose: open=%.3e  gate=%.3e  ratio=%.1fx\n",
              dose_open, dose_gate, dose_open / (dose_gate + 1e-30));

  // Resist blocks: open S/D must receive far more dopant than under the gate.
  CHECK(dose_open > 3.0 * dose_gate);
  // The open region must have actually received dopant.
  CHECK(dose_open > 0);

  std::printf("photo/mask/strip deck tests passed\n");
  return 0;
}
