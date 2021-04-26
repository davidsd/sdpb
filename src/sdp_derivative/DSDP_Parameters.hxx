#pragma once

#include "../Verbosity.hxx"
#include "../sdp_solve.hxx"

#include <El.hpp>
#include <boost/filesystem.hpp>
#include <boost/property_tree/ptree.hpp>

#include <iostream>

struct DSDP_Parameters
{
  size_t procs_per_node, proc_granularity;
  size_t precision;
  Verbosity verbosity;

  boost::filesystem::path primary_sdp_path, secondary_sdp_path, checkpoint_in,
    out_directory, param_path;

  DSDP_Parameters(int argc, char *argv[]);
  bool is_valid() const { return !primary_sdp_path.empty(); }
};

std::ostream &operator<<(std::ostream &os, const DSDP_Parameters &p);

boost::property_tree::ptree to_property_tree(const DSDP_Parameters &p);
