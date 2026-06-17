#pragma once
#include "mesh.hpp"
#include "vec3.hpp"

namespace cp {

// Normalized tetrahedron quality in (0, 1]: q = 6*sqrt(2)*V / l_rms^3, where
// V is the volume and l_rms the root-mean-square edge length. A regular tet
// scores 1; slivers and degenerate cells approach 0. Sign follows the volume.
double tet_quality(const Vec3& a, const Vec3& b, const Vec3& c, const Vec3& d);

struct QualityStats {
  double min_q = 1.0;   // worst cell
  double mean_q = 0.0;  // average over cells
  int n_sliver = 0;     // cells below sliver_thresh
  int worst_cell = -1;
};

QualityStats mesh_quality(const Mesh& m, double sliver_thresh = 0.1);

// Splits edge (a,b) at its midpoint, subdividing every incident tet into two.
// New cells inherit the region of their parent. Returns the new node id, or
// -1 if the edge is not present. The caller must call Mesh::finalize()
// afterward to rebuild faces and geometry.
int edge_split(Mesh& m, int a, int b);

// Laplacian smoothing of interior nodes; boundary and interface nodes are
// pinned. A candidate move is rejected if it would invert any incident tet.
// Performs `iters` sweeps with relaxation `omega`. Returns the number of node
// moves applied. The caller must call Mesh::finalize() afterward.
int laplacian_smooth(Mesh& m, int iters = 5, double omega = 0.5);

}  // namespace cp
