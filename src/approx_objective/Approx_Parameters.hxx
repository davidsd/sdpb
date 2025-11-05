#pragma once

#include "sdp_solve/Solver_Parameters/Memory_Limit.hxx"
#include "sdpb_util/Environment.hxx"
#include "sdpb_util/Verbosity.hxx"

#include <filesystem>
#include <iostream>
#include <boost/property_tree/ptree.hpp>

struct Approx_Parameters
{
  size_t proc_granularity = 1;
  size_t precision;
  Memory_Limit max_memory, max_shared_memory;

  std::filesystem::path sdp_path, new_sdp_path, solution_dir, param_path;

  bool write_solver_state, linear_only;

  Verbosity verbosity;

  // Initialized from Environment - MemAvailable
  String_To_Memory_Limit_Translator memory_limit_translator;

  Approx_Parameters(int argc, char *argv[], const Environment &env);
  [[nodiscard]] bool is_valid() const { return !sdp_path.empty(); }
};

std::ostream &operator<<(std::ostream &os, const Approx_Parameters &p);

boost::property_tree::ptree to_property_tree(const Approx_Parameters &p);
