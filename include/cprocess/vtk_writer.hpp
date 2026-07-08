#pragma once
#include <string>
#include <utility>
#include <vector>

#include "mesh.hpp"

namespace cp {

// Writes the mesh and cell data as a VTK XML unstructured grid (.vtu, ASCII)
// readable by ParaView. Coordinates are written in micrometres for
// convenient viewing; scalar fields are written as-is (cm^-3).
// point_scalars (P3-h): optional node-valued Float64 arrays written as
// <PointData> right after </Cells> and before <CellData>. Empty (default)
// reproduces the original output byte-for-byte.
void write_vtu(
    const std::string& path, const Mesh& mesh,
    const std::vector<std::pair<std::string, const std::vector<double>*>>& scalars,
    const std::vector<std::pair<std::string, const std::vector<int>*>>& int_scalars,
    const std::vector<std::pair<std::string, const std::vector<double>*>>&
        point_scalars = {});

}  // namespace cp
