#pragma once
#include <string>
#include <utility>
#include <vector>

#include "mesh.hpp"

namespace cp {

// Writes the mesh and cell data as a VTK XML unstructured grid (.vtu, ASCII)
// readable by ParaView. Coordinates are written in micrometres for
// convenient viewing; scalar fields are written as-is (cm^-3).
void write_vtu(
    const std::string& path, const Mesh& mesh,
    const std::vector<std::pair<std::string, const std::vector<double>*>>& scalars,
    const std::vector<std::pair<std::string, const std::vector<int>*>>& int_scalars);

}  // namespace cp
