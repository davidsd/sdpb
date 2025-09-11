#pragma once

// #include "SDPA/Polynomial_Vector_Matrix.hxx"

#include <El.hpp>

#include <filesystem>
#include <optional>
#include <vector>

namespace Sdpb::Sdpa
{
  struct SDPA_File_Parse_Result
  {
    // Vector c_1..c_M
    std::vector<El::BigFloat> c_objective;
    // If file is read by several processes,
    // each process saves only some matrices, according to should_parse_block()
    // parsed_blocks is a map: block_index -> blocks of [F_0, F_1,...F_M]
    std::map<size_t, std::vector<El::Matrix<El::BigFloat>>> parsed_blocks;

    SDPA_File_Parse_Result() = default;

    static void validate(const SDPA_File_Parse_Result &result);

    static SDPA_File_Parse_Result
    read(const std::filesystem::path &input_path, bool should_parse_objective,
         const std::function<bool(size_t block_index)> &should_parse_block);

    // Allow moving and prevent accidental copying

    SDPA_File_Parse_Result(const SDPA_File_Parse_Result &other) = delete;
    SDPA_File_Parse_Result(SDPA_File_Parse_Result &&other) noexcept = default;
    SDPA_File_Parse_Result &operator=(const SDPA_File_Parse_Result &other)
      = delete;
    SDPA_File_Parse_Result &operator=(SDPA_File_Parse_Result &&other) noexcept
      = default;
  };
}
