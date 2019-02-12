#pragma once

#include "sdp_convert/Dual_Constraint_Group.hxx"

#include <boost/filesystem.hpp>

#include <vector>

void write_sdpb_input_files(
  const boost::filesystem::path &output_dir, const int &rank,
  const int &num_procs, const std::vector<size_t> &indices,
  const El::BigFloat &objective_const,
  const std::vector<El::BigFloat> &dual_objective_b,
  const std::vector<Dual_Constraint_Group> &dual_constraint_groups);
