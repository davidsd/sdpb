#pragma once

#include "sdpb_util/Verbosity.hxx"
#include "sdp_solve/sdp_solve.hxx"

#include <El.hpp>
#include <filesystem>
#include <boost/property_tree/ptree.hpp>

#include <iostream>

struct SDPB_Parameters
{
  bool no_final_checkpoint;
  size_t proc_granularity;
  bool require_initial_checkpoint = false;
  Write_Solution write_solution;

  Solver_Parameters solver;
  Verbosity verbosity;

  std::filesystem::path sdp_path, out_directory, param_path;

  SDPB_Parameters(int argc, char *argv[]);
  bool is_valid() const { return !sdp_path.empty(); }
};

std::ostream &operator<<(std::ostream &os, const SDPB_Parameters &p);

boost::property_tree::ptree to_property_tree(const SDPB_Parameters &p);
