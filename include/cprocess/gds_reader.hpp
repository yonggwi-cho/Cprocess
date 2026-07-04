#pragma once
#include <string>
#include <utility>
#include <vector>

namespace cp {

// Dependency-free minimal GDSII stream (.gds) binary reader.
//
// Reads BOUNDARY records on the requested layer (layer < 0 means all
// layers) and returns them as closed polygons in centimetres (the GDS
// user coordinates are scaled by the UNITS record's meters-per-DB-unit,
// then converted m -> cm by x100). The trailing duplicate point GDS uses
// to close a polygon is stripped from the returned vertex list.
//
// Hierarchy (SREF/AREF) is NOT expanded -- such records are skipped
// (their payload is read past, not interpreted) along with any other
// unrecognized record type (PATH, TEXT, ...). Malformed files (missing
// HEADER, truncated records, or a nonexistent path) raise
// std::runtime_error.
std::vector<std::vector<std::pair<double, double>>> read_gds(
    const std::string& path, int layer = -1);

}  // namespace cp
