#pragma once

#include "Write_Solution.hxx"
#include "../sdp_solve.hxx"

#include <El.hpp>
#include <boost/filesystem.hpp>
#include <boost/property_tree/ptree.hpp>

#include <iostream>

struct SDPB_Parameters
{
  bool no_final_checkpoint;
  bool require_initial_checkpoint = false;
  size_t precision, procs_per_node, proc_granularity;
  Write_Solution write_solution;

  SDP_Solver_Parameters solver_parameters;
  
  boost::filesystem::path sdp_directory, out_directory, checkpoint_in,
    checkpoint_out, param_file;

  SDPB_Parameters(int argc, char *argv[]);
  bool is_valid() const { return !sdp_directory.empty(); }
};

std::ostream &operator<<(std::ostream &os, const SDPB_Parameters &p);

boost::property_tree::ptree to_property_tree(const SDPB_Parameters &p);
