#pragma once
#include <iosfwd>
#include <map>
#include <string>
#include <vector>

#include "diffusion.hpp"
#include "mc_implant.hpp"
#include "mesh.hpp"

namespace cp {

// Simulation state threaded through a process deck.
struct SimState {
  Mesh mesh;
  bool has_mesh = false;
  // dopant symbol -> per-cell concentration (cm^-3)
  std::map<std::string, std::vector<double>> fields;
  // region tag -> material name ("silicon", "oxide", ...); only silicon
  // cells are implanted/diffused
  std::map<int, std::string> region_material;
  std::vector<DirichletBC> bcs;
  double last_temp = 1273.15;  // K, used for solubility clamping on save

  // ── Physical process stack ────────────────────────────────────────────────
  // Built by 'photo', consumed by 'implant method=mc', destroyed by 'strip'.
  // The stack mesh covers the same (x,y) footprint as `mesh` but extends
  // upward by the resist thickness. MC transport runs in the full stack;
  // the resulting dopant profile is transferred back onto `mesh`.
  bool has_stack = false;
  Mesh stack;                             // substrate + overlayer (Si + resist + vacuum)
  std::vector<int> stack_cell_mat;        // material index per stack cell
  std::vector<TargetMaterial> mat_table;  // [0]=Si, [1]=resist, [2]=vacuum/open
  double stack_resist_z0 = 0;            // z bottom of the resist layer (top of Si)
};

// Runs a process deck (one command per line, '#' comments). Commands:
//   mesh box xmax=1um ymax=1um zmax=2um nx=16 ny=16 nz=32 [xmin= ymin= zmin=]
//   mesh gmsh file=dev.msh [scale=1um]
//   region {tag=N | name=NAME | all} material=silicon|oxide|nitride|poly|gas
//   init species=B conc=1e15 [region=N|NAME]
//   implant species=P dose=1e13 {energy=80keV | rp=0.1um drp=0.04um}
//           [drl=0.03um] [x1= x2= y1= y2=]      (mask window)
//   implant species=B dose=1e13 energy=50keV method=mc
//           [ions=100000] [threads=0] [seed=1] [tilt=7] [rotation=30]
//           [x1= x2= y1= y2=]                   (Monte Carlo / BCA)
//   photo   resist=0.5um [nz=4]      deposit blanket photoresist
//   mask    x1=0.3um x2=0.7um [y1=0 y2=ymax]   expose & develop window
//   strip                             strip all remaining photoresist
//   bc species=P patch=zmax conc=1e20 | bc clear
//   diffuse time=30min temp=1000C [dt=30s] [fieldenh=on|off] [nonortho=on|off]
//   save [file=out.vtu]
//   print
//   stop
void run_deck(std::istream& in, SimState& st, std::ostream& log);

}  // namespace cp
