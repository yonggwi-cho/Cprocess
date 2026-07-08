#pragma once
#include <string>
#include <vector>

#include "deck.hpp"
#include "mesh.hpp"

namespace cp {

// Volume-weighted cell->node averaging (P3-h):
//   val(n) = sum_{c: n in mesh.cells[c], exclude[c]==0} cell_vol[c]*f[c]
//          / sum_{c: n in mesh.cells[c], exclude[c]==0} cell_vol[c]
// exclude[c] != 0 marks cells skipped (e.g. gas); a node adjacent only to
// excluded cells (or to no cells at all) gets 0.0.
std::vector<double> cell_to_node(const Mesh& mesh,
                                  const std::vector<double>& f,
                                  const std::vector<char>& exclude);

// Writes <prefix>.vtu (cell data identical to proc::save plus node-averaged
// dopant/active/NetDoping point data) and <prefix>.meta.json (region/
// material table, boundary patch names, unit system, species list, ND-NA
// sign convention). Throws std::runtime_error on I/O failure.
void write_device(const SimState& st, const std::string& prefix);

}  // namespace cp
