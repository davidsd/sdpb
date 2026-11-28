#pragma once

#include "SDPA_Block_Structure.hxx"
#include "SDPA_File_Parse_Result.hxx"

#include <filesystem>
#include <functional>

namespace Sdpb::Sdpa
{
SDPA_Block_Structure
read_block_structure(const std::filesystem::path &input_file);

SDPA_File_Parse_Result read_sdpa_objective(const std::filesystem::path &input_file);

// Visit the lower triangle without materializing dense F matrices.
// Indices are zero-based, with diagonal blocks split into scalar blocks.
// Sparse files are cut into slices at line ends, and the slices are parsed on
// OMP_NUM_THREADS threads: visit runs concurrently for different slots
// (0 <= slot < sdpa_visit_slots()), never concurrently for the same slot.
// order increases along the file, so sorting by it recovers file order (later
// lines overwrite earlier ones).
using SDPA_Entry_Visitor = std::function<void(
  size_t slot, uint64_t order, size_t block, size_t matrix, size_t i, size_t j,
  const El::BigFloat &value)>;
size_t sdpa_visit_slots();
void visit_sdpa_lower(const std::filesystem::path &input_file,
                      const SDPA_Entry_Visitor &visit);

SDPA_File_Parse_Result
read_sdpa(const std::filesystem::path &input_file,
          const std::function<bool(size_t matrix_index)> &should_parse_block);
}