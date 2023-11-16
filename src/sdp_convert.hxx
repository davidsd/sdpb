#pragma once

#include "Boost_Float.hxx"
#include "sdp_convert/Block_File_Format.hxx"
#include "sdp_convert/Dual_Constraint_Group.hxx"
#include "sdp_convert/write_vector.hxx"
#include "Timers.hxx"

#include <filesystem>

#include <vector>

void write_sdpb_input_files(
  const std::filesystem::path &output_path, Block_File_Format output_format,
  const int &rank, const size_t &num_blocks,
  const std::vector<std::string> &command_arguments,
  const El::BigFloat &objective_const,
  const std::vector<El::BigFloat> &dual_objective_b,
  const std::vector<Dual_Constraint_Group> &dual_constraint_groups,
  Timers &timers, const bool debug);
