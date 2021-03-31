#pragma once

#include "../Verbosity.hxx"
#include "../sdp_solve.hxx"

#include <El.hpp>
#include <boost/filesystem.hpp>
#include <boost/property_tree/ptree.hpp>

#include <iostream>

struct Outer_Parameters
{
  bool no_final_checkpoint;
  bool require_initial_checkpoint = false;
  Write_Solution write_solution;

  Solver_Parameters solver;
  Verbosity verbosity;

  boost::filesystem::path sdp_path, points_path, out_directory, param_path;

  Outer_Parameters(int argc, char *argv[]);
  bool is_valid() const { return !sdp_path.empty(); }
};

std::ostream &operator<<(std::ostream &os, const Outer_Parameters &p);

boost::property_tree::ptree to_property_tree(const Outer_Parameters &p);
