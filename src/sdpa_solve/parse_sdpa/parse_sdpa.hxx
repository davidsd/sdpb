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

// Visit the lower triangle in file order without materializing dense F matrices.
// Indices are zero-based, with diagonal blocks split into scalar blocks.
using SDPA_Entry_Visitor = std::function<void(size_t block, size_t matrix,
                                               size_t i, size_t j,
                                               const El::BigFloat &value)>;
void visit_sdpa_lower(const std::filesystem::path &input_file,
                      const SDPA_Entry_Visitor &visit);

SDPA_File_Parse_Result
read_sdpa(const std::filesystem::path &input_file,
          const std::function<bool(size_t matrix_index)> &should_parse_block);
}