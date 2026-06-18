// examples/nmos_locos_flow.cpp
//
// Integrated example: simplified LOCOS-style NMOS process flow
//
// Process steps (all units match Cprocess conventions):
//   1. Initial Si substrate: 1 um x 1 um x 0.5 um, p-type (boron background)
//   2. Shallow boron channel implant (p-well punch-through stop)
//   3. Phosphorus source/drain implant with photoresist mask
//   4. Rapid thermal anneal: diffusion at 1000°C × 30 s
//   5. Gate oxide simulation: Deal-Grove dry oxidation at 950°C × 10 min
//   6. ALE mesh motion: grow SiO2 layer on top surface
//   7. Level-set region tagging: Si vs SiO2 via φ sign
//   8. Mesh quality report + Laplacian smoothing
//   9. Viscoelastic stress in the oxide (Maxwell relaxation)
//  10. Field transfer: dopant profile onto a refined mesh
//  11. VTK output: boron/phosphorus concentrations, stress, level-set
//
// Build (from Cprocess root):
//   cmake -B build -DCMAKE_BUILD_TYPE=Release
//   cmake --build build -j
//   ./build/nmos_locos_flow
//
// Output: nmos_locos_*.vtu  (open with ParaView)

#include <array>
#include <cmath>
#include <cstdio>
#include <vector>
#include <utility>

#include "cprocess/ale_mover.hpp"
#include "cprocess/diffusion.hpp"
#include "cprocess/field_transfer.hpp"
#include "cprocess/levelset.hpp"
#include "cprocess/materials.hpp"
#include "cprocess/mc_implant.hpp"
#include "cprocess/mechanics.hpp"
#include "cprocess/mesh.hpp"
#include "cprocess/oxidation.hpp"
#include "cprocess/remesh.hpp"
#include "cprocess/topology.hpp"
#include "cprocess/vtk_writer.hpp"

using namespace cp;

// ── helpers ──────────────────────────────────────────────────────────────────

static void print_step(int n, const char* label) {
  std::printf("\n══ Step %d: %s ══\n", n, label);
}

static void print_stats(const char* name, const std::vector<double>& f,
                        const std::vector<double>& vol) {
  double peak = 0, total = 0;
  for (std::size_t i = 0; i < f.size(); ++i) {
    if (f[i] > peak) peak = f[i];
    total += f[i] * vol[i];
  }
  std::printf("  %-10s  peak=%8.3e cm⁻³  integral=%8.3e cm⁻²\n",
              name, peak, total / 1e-4 /* z-thickness 0.5 um in cm */);
}

// ── main ─────────────────────────────────────────────────────────────────────

int main() {
  std::printf("╔══════════════════════════════════════════════════════╗\n");
  std::printf("║   Cprocess  –  NMOS LOCOS Integration Example       ║\n");
  std::printf("╚══════════════════════════════════════════════════════╝\n");

  // ── 1. Substrate mesh ─────────────────────────────────────────────────────
  print_step(1, "Substrate mesh");
  // 1 um × 1 um × 0.5 um Si, coarse 8×8×4 hex → 6×8×8×4 = 1536 tets.
  // All lengths in cm; 1 um = 1e-4 cm.
  const double Lx = 1e-4, Ly = 1e-4, Lz = 0.5e-4;
  const int NX = 8, NY = 8, NZ = 4;
  Mesh mesh = make_box_mesh(0, Lx, 0, Ly, 0, Lz, NX, NY, NZ);
  const int nc = static_cast<int>(mesh.cells.size());
  const int nn = static_cast<int>(mesh.nodes.size());
  std::printf("  mesh: %d tets, %d nodes\n", nc, nn);

  // Region 1 = Si substrate.
  for (int& r : mesh.cell_region) r = 1;
  mesh.region_names[1] = "Si";
  mesh.region_names[2] = "SiO2";

  // Per-cell concentration fields (cm⁻³).
  std::vector<double> c_boron(nc, 0.0);
  std::vector<double> c_phos(nc, 0.0);

  // silicon mask: all cells are Si at the start.
  std::vector<char> si_mask(nc, 1);

  // ── 2. Channel boron implant (background p-well) ──────────────────────────
  print_step(2, "Boron channel implant");
  {
    const Dopant* B = find_dopant("B");
    if (!B) { std::fprintf(stderr, "ERROR: boron not found\n"); return 1; }

    McImplantParams p;
    p.dopant      = B;
    p.dose        = 5e12;        // cm⁻²
    p.energy_kev  = 40;
    p.tilt_deg    = 7;
    p.ions        = 50000;
    p.channeling  = true;
    p.seed        = 42;

    const auto stats = apply_mc_implant(mesh, si_mask, p, c_boron, nullptr);
    std::printf("  B 40keV: Rp=%.1f nm  deposited=%lld  backscattered=%lld\n",
                stats.rp * 1e7, stats.deposited, stats.backscattered);
    print_stats("boron", c_boron, mesh.cell_vol);
  }

  // ── 3. Phosphorus S/D implant through a PHYSICAL patterned resist ─────────
  print_step(3, "Phosphorus S/D implant (physical photoresist gate)");
  {
    const Dopant* P = find_dopant("P");
    if (!P) { std::fprintf(stderr, "ERROR: phosphorus not found\n"); return 1; }

    // Build a process stack on top of the substrate: Si [0,0.5 um] + a 0.5 um
    // tall overlayer. A photoresist gate covers the channel x ∈ [0.35,0.65 um];
    // the developed source/drain openings are filled with near-vacuum so ions
    // reach the real silicon surface. The whole top is irradiated (no window) —
    // masking is entirely physical: ions stop inside the resist or straggle
    // laterally under its edge.
    const double t_over = 0.5e-4;            // overlayer thickness
    const double Lz_stack = Lz + t_over;
    const double gate_x1 = 0.35e-4, gate_x2 = 0.65e-4;
    Mesh stack = make_box_mesh(0, Lx, 0, Ly, 0, Lz_stack, NX, NY, NZ + 4);
    const int nc_s = static_cast<int>(stack.cells.size());

    // Material map: 0 = Si, 1 = photoresist, 2 = vacuum (open S/D).
    enum { MAT_SI = 0, MAT_RESIST = 1, MAT_VAC = 2 };
    std::vector<int> mat(nc_s, MAT_SI);
    std::vector<char> stack_si(nc_s, 1);
    for (int ci = 0; ci < nc_s; ++ci) {
      const Vec3& c = stack.cell_cent[ci];
      if (c.z <= Lz) continue;             // silicon substrate
      stack_si[ci] = 0;                    // overlayer never accumulates dopant
      const bool under_gate = (c.x >= gate_x1 && c.x <= gate_x2);
      mat[ci] = under_gate ? MAT_RESIST : MAT_VAC;
    }

    McImplantParams p;
    p.dopant         = P;
    p.dose           = 5e15;     // cm⁻²
    p.energy_kev     = 30;
    p.tilt_deg       = 0;
    p.ions           = 80000;
    p.channeling     = true;
    p.seed           = 7;
    p.material_table = {target_silicon(), target_photoresist(), target_vacuum()};
    p.cell_material  = &mat;     // full-surface irradiation, physical mask

    std::vector<double> c_phos_stack(nc_s, 0.0);
    const auto stats = apply_mc_implant(stack, stack_si, p, c_phos_stack, nullptr);
    std::printf("  P 30keV: deposited(Si)=%lld  stopped_in_overlayer=%lld\n",
                stats.deposited, stats.in_mask);

    // Lateral check: dopant under the gate vs in the open S/D, in silicon.
    double q_gate = 0, q_open = 0;
    for (int ci = 0; ci < nc_s; ++ci) {
      if (!stack_si[ci]) continue;
      const Vec3& c = stack.cell_cent[ci];
      const double q = c_phos_stack[ci] * stack.cell_vol[ci];
      if (c.x >= gate_x1 && c.x <= gate_x2) q_gate += q; else q_open += q;
    }
    std::printf("  silicon P dose: open=%.3e  under-gate=%.3e  blocking=%.0fx\n",
                q_open, q_gate, q_open / (q_gate + 1e-30));

    // Transfer the phosphorus from the stack's silicon onto the working mesh
    // (P4 conservative nearest-cell transfer).
    c_phos = transfer_field_nearest(stack, c_phos_stack, mesh);
    print_stats("phosphorus", c_phos, mesh.cell_vol);
  }

  // ── 4. RTA diffusion: 1000°C × 30 s ──────────────────────────────────────
  print_step(4, "Diffusion anneal 1000°C 30 s");
  {
    const Dopant* B = find_dopant("B");
    const Dopant* P = find_dopant("P");

    DiffuseOpts opts;
    opts.temp      = 1273.15;   // K
    opts.time      = 30.0;      // s
    opts.dt        = 0.5;       // s per step
    opts.field_enh = true;
    opts.verbosity = 0;

    DiffusionSolver solver(mesh, si_mask, nullptr);

    std::vector<SpeciesField> fields = {
      {B, &c_boron},
      {P, &c_phos},
    };
    solver.run(fields, {}, opts);

    std::printf("  after anneal:\n");
    print_stats("boron", c_boron, mesh.cell_vol);
    print_stats("phosphorus", c_phos, mesh.cell_vol);
  }

  // Write checkpoint after implant+diffuse.
  {
    std::vector<std::pair<std::string, const std::vector<double>*>> scalars = {
      {"boron",     &c_boron},
      {"phosphorus",&c_phos},
    };
    write_vtu("nmos_locos_01_after_diffuse.vtu", mesh, scalars, {});
    std::printf("  wrote nmos_locos_01_after_diffuse.vtu\n");
  }

  // ── 5. Gate oxide thickness (Deal-Grove) ──────────────────────────────────
  print_step(5, "Gate oxide: Deal-Grove dry 950°C 10 min");
  const double tox_um = deal_grove_step(0.0, 10.0, 950.0, false);
  std::printf("  tox = %.2f nm (%.4f um)\n", tox_um * 1000.0, tox_um);
  // SiO2 grows by 2.2× the consumed Si. Convert to cm for mesh motion.
  const double tox_cm = tox_um * 1e-4;

  // ── 6. ALE: grow top surface by tox_cm ───────────────────────────────────
  print_step(6, "ALE mesh motion: oxide growth on zmax");
  {
    MeshTopology topo;
    topo.build(mesh);
    const auto normals = compute_node_normals(mesh, topo);

    // Save old volumes for field rescaling.
    const std::vector<double> old_vol = mesh.cell_vol;

    // Only move zmax boundary nodes (top surface).
    const double zmax_before = mesh.bbox().hi.z;
    auto mask = [&](int i) -> bool {
      return !topo.node_boundary[i] || (mesh.nodes[i].z < zmax_before * 0.999);
    };

    // v_n = tox_cm, dt = 1 (single shot: move by tox_cm).
    const auto res = ale_move(mesh, topo, normals, tox_cm, 1.0, mask, 3);
    mesh.finalize();

    const double zmax_after = mesh.bbox().hi.z;
    std::printf("  n_moved=%d  zmax %.4f → %.4f um\n",
                res.n_moved, zmax_before * 1e4, zmax_after * 1e4);

    // Rescale dopant concentrations to conserve mass.
    std::vector<std::vector<double>> fields = {c_boron, c_phos};
    rescale_fields_for_volume_change(fields, old_vol, mesh.cell_vol);
    c_boron = fields[0];
    c_phos  = fields[1];
  }

  // ── 7. Level-set: tag SiO2 region ────────────────────────────────────────
  print_step(7, "Level-set: tag oxide cells as SiO2");
  {
    // After ALE, cells whose nodes were displaced (max node z > Lz) represent
    // the stretched top layer = SiO2. Tag via max node z > Lz - epsilon.
    const double z_thresh = Lz - Lz / (2.0 * NZ);  // halfway into top layer
    for (int ci = 0; ci < nc; ++ci) {
      double zmax_cell = 0;
      for (int v : mesh.cells[ci]) zmax_cell = std::max(zmax_cell, mesh.nodes[v].z);
      mesh.cell_region[ci] = (zmax_cell > z_thresh) ? 2 : 1;
    }

    // Initialize φ: Si = inside (region 1), SiO2 = outside.
    const auto phi = levelset_init(mesh, 1);

    // Update si_mask: only region-1 cells diffuse.
    for (int ci = 0; ci < nc; ++ci)
      si_mask[ci] = (mesh.cell_region[ci] == 1) ? 1 : 0;

    int n_si = 0, n_ox = 0;
    for (int r : mesh.cell_region) { if (r == 1) ++n_si; else ++n_ox; }
    std::printf("  Si cells=%d  SiO2 cells=%d\n", n_si, n_ox);

    // Write level-set on nodes for visualization.
    // (Store on cells as average of node values for VTK.)
    std::vector<double> phi_cell(nc);
    for (int ci = 0; ci < nc; ++ci) {
      for (int v : mesh.cells[ci]) phi_cell[ci] += phi[v];
      phi_cell[ci] /= 4.0;
    }

    std::vector<std::pair<std::string, const std::vector<double>*>> scalars = {
      {"boron",     &c_boron},
      {"phosphorus",&c_phos},
      {"levelset",  &phi_cell},
    };
    std::vector<std::pair<std::string, const std::vector<int>*>> int_s = {
      {"region", &mesh.cell_region},
    };
    write_vtu("nmos_locos_02_after_ale.vtu", mesh, scalars, int_s);
    std::printf("  wrote nmos_locos_02_after_ale.vtu\n");
  }

  // ── 8. Mesh quality ───────────────────────────────────────────────────────
  print_step(8, "Mesh quality + Laplacian smoothing");
  {
    QualityStats qs0 = mesh_quality(mesh);
    std::printf("  before smooth: min_q=%.4f  mean_q=%.4f  slivers=%d\n",
                qs0.min_q, qs0.mean_q, qs0.n_sliver);

    const int moved = laplacian_smooth(mesh, 5, 0.5);
    mesh.finalize();

    QualityStats qs1 = mesh_quality(mesh);
    std::printf("  after  smooth: min_q=%.4f  mean_q=%.4f  moved=%d nodes\n",
                qs1.min_q, qs1.mean_q, moved);
  }

  // ── 9. Viscoelastic stress in SiO2 ───────────────────────────────────────
  print_step(9, "Viscoelastic oxide stress (Maxwell)");
  {
    // SiO2 parameters: E=70 GPa, nu=0.17, tau ~ 100 min at 950°C.
    ViscoElasticParams sio2;
    sio2.E_GPa       = 70.0;
    sio2.nu          = 0.17;
    sio2.tau_relax_min = 100.0;

    std::vector<StressState> stress(nc);
    // Thermal mismatch strain: SiO2 CTE 0.5e-6/K, Si 2.6e-6/K, ΔT = 725 K
    const double d_alpha = (0.5e-6 - 2.6e-6);
    const double dT = 725.0;
    const double eps_th = d_alpha * dT;  // biaxial mismatch in x and y

    std::vector<double> dstrain(nc * 6, 0.0);
    for (int ci = 0; ci < nc; ++ci) {
      if (mesh.cell_region[ci] != 2) continue;
      dstrain[ci * 6 + 0] = eps_th;  // εxx
      dstrain[ci * 6 + 1] = eps_th;  // εyy
    }

    // Single time step: anneal duration (10 min → already happened; model relaxation).
    mechanics_step(stress, dstrain, 10.0, sio2);

    // Collect σxx for output.
    std::vector<double> sxx(nc, 0.0);
    double peak_stress_MPa = 0;
    for (int ci = 0; ci < nc; ++ci) {
      sxx[ci] = stress[ci].s[0];
      if (std::fabs(sxx[ci]) > std::fabs(peak_stress_MPa))
        peak_stress_MPa = sxx[ci];
    }
    std::printf("  peak σxx = %.1f MPa (in SiO2)\n", peak_stress_MPa);

    std::vector<std::pair<std::string, const std::vector<double>*>> scalars = {
      {"boron",     &c_boron},
      {"phosphorus",&c_phos},
      {"stress_xx", &sxx},
    };
    std::vector<std::pair<std::string, const std::vector<int>*>> int_s = {
      {"region", &mesh.cell_region},
    };
    write_vtu("nmos_locos_03_stress.vtu", mesh, scalars, int_s);
    std::printf("  wrote nmos_locos_03_stress.vtu\n");
  }

  // ── 10. Field transfer to refined mesh ───────────────────────────────────
  print_step(10, "Field transfer → refined mesh");
  {
    const double new_Lz = mesh.bbox().hi.z;
    Mesh mesh_fine = make_box_mesh(0, Lx, 0, Ly, 0, new_Lz, 12, 12, 6);
    const int nc_fine = static_cast<int>(mesh_fine.cells.size());

    const auto c_boron_fine = transfer_field_nearest(mesh, c_boron, mesh_fine);
    const auto c_phos_fine  = transfer_field_nearest(mesh, c_phos,  mesh_fine);

    const double err_B = check_mass_conservation(mesh, c_boron, mesh_fine, c_boron_fine);
    const double err_P = check_mass_conservation(mesh, c_phos,  mesh_fine, c_phos_fine);
    std::printf("  fine mesh: %d tets\n", nc_fine);
    std::printf("  mass conservation error: B=%.2f%%  P=%.2f%%\n",
                err_B * 100.0, err_P * 100.0);

    std::vector<std::pair<std::string, const std::vector<double>*>> scalars = {
      {"boron",     &c_boron_fine},
      {"phosphorus",&c_phos_fine},
    };
    write_vtu("nmos_locos_04_fine.vtu", mesh_fine, scalars, {});
    std::printf("  wrote nmos_locos_04_fine.vtu\n");
  }

  // ── Summary ───────────────────────────────────────────────────────────────
  std::printf("\n╔══════════════════════════════════════════════════════╗\n");
  std::printf("║  All steps completed successfully.                   ║\n");
  std::printf("║  Output files:                                       ║\n");
  std::printf("║    nmos_locos_01_after_diffuse.vtu                   ║\n");
  std::printf("║    nmos_locos_02_after_ale.vtu                       ║\n");
  std::printf("║    nmos_locos_03_stress.vtu                          ║\n");
  std::printf("║    nmos_locos_04_fine.vtu                            ║\n");
  std::printf("║  Open with ParaView to inspect 3D profiles.          ║\n");
  std::printf("╚══════════════════════════════════════════════════════╝\n");
  return 0;
}
