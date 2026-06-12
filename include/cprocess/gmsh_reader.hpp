#pragma once
#include <iosfwd>
#include <string>

#include "mesh.hpp"

namespace cp {

// Reads a Gmsh ASCII .msh file (format 2.2 or 4.1), first-order tets only.
//
//  - physical volumes  -> cell regions (tag + name)
//  - physical surfaces -> boundary patches (tagged triangles)
//  - untagged boundary faces are collected into a patch named "default"
//  - node coordinates are multiplied by `scale` (mesh unit in cm,
//    e.g. 1e-4 for a mesh drawn in micrometres)
Mesh read_gmsh(const std::string& path, double scale, std::ostream* log = nullptr);

}  // namespace cp
