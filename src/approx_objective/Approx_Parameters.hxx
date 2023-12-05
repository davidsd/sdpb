#pragma once

#include <El.hpp>
#include <filesystem>
#include <iostream>
#include <boost/property_tree/ptree.hpp>

struct Approx_Parameters
{
  size_t procs_per_node, proc_granularity;
  size_t precision;

  std::filesystem::path sdp_path, new_sdp_path, solution_dir, param_path;

  bool write_solver_state, linear_only;

  Approx_Parameters(int argc, char *argv[]);
  bool is_valid() const { return !sdp_path.empty(); }
};

std::ostream &operator<<(std::ostream &os, const Approx_Parameters &p);

boost::property_tree::ptree to_property_tree(const Approx_Parameters &p);
