#pragma once

#include "../Verbosity.hxx"
#include "../sdp_solve.hxx"

#include <El.hpp>
#include <boost/filesystem.hpp>
#include <boost/property_tree/ptree.hpp>

#include <iostream>

struct Outer_Parameters
{
  bool require_initial_checkpoint = false,
    use_svd = false;
  Write_Solution write_solution;

  El::BigFloat duality_gap_reduction, mesh_threshold;
  Solver_Parameters solver;
  Verbosity verbosity;

  boost::filesystem::path functions_path, points_path, output_path, param_path;

  Outer_Parameters(int argc, char *argv[]);
  bool is_valid() const { return !functions_path.empty(); }
};

std::ostream &operator<<(std::ostream &os, const Outer_Parameters &p);

boost::property_tree::ptree to_property_tree(const Outer_Parameters &p);
