#pragma once

#include "../Block_File_Format.hxx"
#include "sdpb_util/Verbosity.hxx"

#include <boost/program_options.hpp>
#include <boost/property_tree/ptree.hpp>

#include <filesystem>

struct Pmp2sdp_Parameters
{
  int precision;
  std::filesystem::path input_file;
  std::filesystem::path output_path;
  Block_File_Format output_format;
  bool zip = false;
  Verbosity verbosity;

  std::vector<std::string> command_arguments;


  Pmp2sdp_Parameters(int argc, char **argv);
  [[nodiscard]] bool is_valid() const;
};

std::ostream &operator<<(std::ostream &os, const Pmp2sdp_Parameters &p);

boost::property_tree::ptree to_property_tree(const Pmp2sdp_Parameters &p);