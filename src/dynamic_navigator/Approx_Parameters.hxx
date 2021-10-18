#pragma once

#include "../sdp_solve.hxx"

#include <El.hpp>
#include <boost/filesystem.hpp>
#include <boost/property_tree/ptree.hpp>

#include <iostream>

struct Approx_Parameters
{
  size_t procs_per_node, proc_granularity;
  size_t precision;

  boost::filesystem::path sdp_path, new_sdp_path, solution_dir,
    param_path;

  bool write_solver_state, linear_only;

  Approx_Parameters(int argc, char *argv[]);
  bool is_valid() const { return !sdp_path.empty(); }
};

std::ostream &operator<<(std::ostream &os, const Approx_Parameters &p);

boost::property_tree::ptree to_property_tree(const Approx_Parameters &p);
