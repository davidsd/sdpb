#pragma once

#include "sdpb_util/Verbosity.hxx"
#include "dynamical_solve/dynamical_solve.hxx"

#include <filesystem>
#include <boost/property_tree/ptree.hpp>

struct Dynamical_Parameters
{
  bool no_final_checkpoint;
  size_t procs_per_node, proc_granularity;
  bool require_initial_checkpoint = false;
  Write_Solution write_solution;

  Dynamical_Solver_Parameters solver;
  Verbosity verbosity;

  std::filesystem::path sdp_path, out_directory, param_path;

  Dynamical_Parameters(int argc, char *argv[]);
  bool is_valid() const { return !sdp_path.empty(); }
};

std::ostream &operator<<(std::ostream &os, const Dynamical_Parameters &p);

boost::property_tree::ptree to_property_tree(const Dynamical_Parameters &p);
