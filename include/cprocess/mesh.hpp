#pragma once
#include <array>
#include <cstddef>
#include <map>
#include <string>
#include <unordered_map>
#include <vector>

#include "vec3.hpp"

namespace cp {

// Cell-centered unstructured tetrahedral mesh.
//
// Conventions (all lengths in cm):
//  - cells are tetrahedra with positive signed volume (finalize() reorients)
//  - a face is a triangle shared by at most two cells; its area vector S
//    points from `owner` to `neigh`, outward for boundary faces
//  - boundary faces carry a patch index into patch_names; internal faces
//    have patch == -1
struct Face {
  std::array<int, 3> n{};  // node ids, oriented so S points owner -> neigh
  int owner = -1;
  int neigh = -1;          // -1 for boundary faces
  int patch = -1;          // boundary patch index, -1 for internal faces
  Vec3 S;                  // area vector
  Vec3 c;                  // centroid
};

struct TriKeyHash {
  std::size_t operator()(const std::array<int, 3>& a) const {
    std::size_t h = 1469598103934665603ull;
    for (int v : a) {
      h ^= static_cast<std::size_t>(v);
      h *= 1099511628211ull;
    }
    return h;
  }
};

struct BBox {
  Vec3 lo, hi;
};

struct Mesh {
  std::vector<Vec3> nodes;
  std::vector<std::array<int, 4>> cells;
  std::vector<int> cell_region;               // region tag per cell
  std::map<int, std::string> region_names;    // region tag -> name
  std::vector<std::string> patch_names;       // boundary patch names
  std::vector<Face> faces;
  std::vector<double> cell_vol;
  std::vector<Vec3> cell_cent;
  // sorted node triple -> face index, built by finalize()
  std::unordered_map<std::array<int, 3>, int, TriKeyHash> face_lookup;

  // Orients cells to positive volume, builds the face list and all
  // geometric quantities. Must be called once after nodes/cells are set.
  void finalize();

  int add_patch(const std::string& name);            // returns existing id if present
  int find_patch(const std::string& name) const;     // -1 if absent
  BBox bbox() const;
  double total_volume() const;
  // min over internal faces of cos(angle between S and the P->N vector);
  // must be > 0 for a monotone two-point flux scheme
  double min_orthogonality() const;
  std::vector<int> region_tags() const;              // distinct, sorted
};

// Axis-aligned box meshed with a conforming 6-tet (Kuhn) subdivision of an
// nx*ny*nz hex grid. Boundary patches: xmin,xmax,ymin,ymax,zmin,zmax.
Mesh make_box_mesh(double x0, double x1, double y0, double y1, double z0,
                   double z1, int nx, int ny, int nz);

}  // namespace cp
