#pragma once

#include "sdpb_util/Verbosity.hxx"
#include "sdp_solve/sdp_solve.hxx"

#include <El.hpp>
#include <filesystem>
#include <boost/property_tree/ptree.hpp>

#include <iostream>
#include <string>
#include <vector>

struct SDPB_Parameters
{
  bool no_final_checkpoint;
  size_t proc_granularity = 1;
  bool require_initial_checkpoint = false;
  Write_Solution write_solution;

  Solver_Parameters solver;
  Verbosity verbosity;

  std::filesystem::path sdp_path, out_directory, param_path;

  // Resident mode: read commands from this directory instead of solving once
  // (see src/sdpb/resident/Readme.md). Empty: the usual single solve.
  std::filesystem::path control_dir;
  // Resident mode: exit after this many seconds without a command (0: never)
  int64_t control_idle_timeout = 0;

  SDPB_Parameters(int argc, char *argv[], const Environment &env);
  // Parse several argument lists (without the program name). An option is
  // taken from the first list that has it. With throw_on_error, errors throw
  // instead of aborting MPI (used by the resident mode for its solve
  // commands).
  SDPB_Parameters(const std::vector<std::vector<std::string>> &arg_lists,
                  bool throw_on_error, const Environment &env);
  bool is_valid() const { return !sdp_path.empty(); }
  // .dat, .dat-s, .dat.gz or .dat-s.gz
  bool is_sdpa() const;
};

std::ostream &operator<<(std::ostream &os, const SDPB_Parameters &p);

boost::property_tree::ptree to_property_tree(const SDPB_Parameters &p);
