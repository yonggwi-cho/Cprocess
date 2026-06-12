#pragma once
#include <iosfwd>
#include <map>
#include <string>
#include <vector>

#include "diffusion.hpp"
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
};

// Runs a process deck (one command per line, '#' comments). Commands:
//   mesh box xmax=1um ymax=1um zmax=2um nx=16 ny=16 nz=32 [xmin= ymin= zmin=]
//   mesh gmsh file=dev.msh [scale=1um]
//   region {tag=N | name=NAME | all} material=silicon|oxide|nitride|poly|gas
//   init species=B conc=1e15 [region=N|NAME]
//   implant species=P dose=1e13 {energy=80keV | rp=0.1um drp=0.04um}
//           [drl=0.03um] [x1= x2= y1= y2=]      (mask window)
//   bc species=P patch=zmax conc=1e20 | bc clear
//   diffuse time=30min temp=1000C [dt=30s] [fieldenh=on|off] [nonortho=on|off]
//   save [file=out.vtu]
//   print
//   stop
void run_deck(std::istream& in, SimState& st, std::ostream& log);

}  // namespace cp
