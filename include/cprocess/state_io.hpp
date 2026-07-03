#pragma once
#include <string>

#include "deck.hpp"

namespace cp {

// CPRC1 v1 binary state format (little-endian, host representation, x86/ARM64
// only — no endian conversion). Layout, in order:
//   magic "CPRC1" (char[5]), version int64=1,
//   n_nodes int64, node xyz double x 3*n_nodes,
//   n_cells int64, cell node ids int64 x 4*n_cells,
//   cell_region int64 x n_cells,
//   n_region_names, then (tag int64, name_len int64, name bytes) each,
//   n_patch_names, then (name_len int64, name bytes) each,
//   n_region_material, then (tag int64, mat_len int64, mat bytes) each,
//   n_fields, then (sym_len int64, sym bytes, values double x n_cells) each,
//   n_bcs, then (sym_len int64, sym bytes, patch int64, conc double) each,
//   last_temp double.
// The photoresist stack (has_stack / stack / stack_cell_mat / mat_table /
// stack_resist_z0) is never serialized: save_state() throws if has_stack is
// set (strip it first). layer_stack (added to SimState after this format was
// designed) is likewise NOT part of the v1 format; a future format bump that
// adds it should become version 2.
void write_state(const SimState& st, const std::string& path);
void read_state(SimState& st, const std::string& path);

}  // namespace cp
