#pragma once

#include "pmp/Polynomial_Matrix_Program.hxx"
#include "../Dual_Constraint_Group.hxx"
#include "sdpb_util/Timers/Timers.hxx"

#include <El.hpp>

#include <vector>

#include <boost/noncopyable.hpp>

// Data to be written to sdp directory
// see eq. (2.2) in SDPB Manual
struct Output_SDP : boost::noncopyable
{
  // objective_const and dual_objective_b are copied across all processes
  // b_0
  El::BigFloat objective_const;
  // b_1..b_N
  std::vector<El::BigFloat> dual_objective_b;
  // n_0..n_N
  std::vector<El::BigFloat> normalization;
  // Total number of blocks for all processes
  size_t num_blocks = 0;
  // Each dual constraint group is written to  block_info_XXX.json and block_data_XXX.bin.
  // NB: if there are several processes, each of them owns only some of the blocks.
  std::vector<Dual_Constraint_Group> dual_constraint_groups;
  // Command-line arguments
  std::vector<std::string> command_arguments;

  Output_SDP(const Polynomial_Matrix_Program &pmp, const std::vector<std::string> &command_arguments, Timers& timers);
};
