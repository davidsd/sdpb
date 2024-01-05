#pragma once

#include "sdpb_util/Boost_Float.hxx"
#include "Block_File_Format.hxx"
#include "Dual_Constraint_Group.hxx"
#include "write_vector.hxx"
#include "sdpb_util/Timers/Timers.hxx"

#include <filesystem>

#include <vector>

void write_sdpb_input_files(
  const std::filesystem::path &output_path, Block_File_Format output_format,
  const size_t &num_blocks, const std::vector<std::string> &command_arguments,
  const El::BigFloat &objective_const,
  const std::vector<El::BigFloat> &dual_objective_b,
  const std::vector<Dual_Constraint_Group> &dual_constraint_groups,
  Timers &timers, const bool debug);
