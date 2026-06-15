#include "cprocess/implant.hpp"

#include <cmath>
#include <stdexcept>

namespace cp {

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

  conc.resize(mesh.cells.size(), 0.0);
  double atoms = 0;
  for (std::size_t ci = 0; ci < mesh.cells.size(); ++ci) {
    if (!mask.empty() && !mask[ci]) continue;
    const Vec3& c = mesh.cell_cent[ci];
    const double d = ztop - c.z;
    double v = peak * std::exp(-((d - p.rp) * (d - p.rp)) / (s2v * s2v));
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
