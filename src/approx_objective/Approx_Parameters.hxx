#pragma once

#include "sdpb_util/Verbosity.hxx"

#include <filesystem>
#include <iostream>
#include <boost/property_tree/ptree.hpp>

struct Approx_Parameters
{
  size_t proc_granularity = 1;
  size_t precision;
  size_t max_shared_memory_bytes;

  std::filesystem::path sdp_path, new_sdp_path, solution_dir, param_path;

  bool write_solver_state, linear_only;

  Verbosity verbosity;

  Approx_Parameters(int argc, char *argv[]);
  [[nodiscard]] bool is_valid() const { return !sdp_path.empty(); }
};

std::ostream &operator<<(std::ostream &os, const Approx_Parameters &p);

boost::property_tree::ptree to_property_tree(const Approx_Parameters &p);
