#pragma once

#include "SDPA_Block_Structure.hxx"
#include "SDPA_File_Parse_Result.hxx"

#include <filesystem>
#include <functional>

namespace Sdpb::Sdpa
{
SDPA_Block_Structure
read_block_structure(const std::filesystem::path &input_file);

SDPA_File_Parse_Result
read_sdpa(const std::filesystem::path &input_file,
          const std::function<bool(size_t matrix_index)> &should_parse_block);
}