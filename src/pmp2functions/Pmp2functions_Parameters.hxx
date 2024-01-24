#pragma once

#include "sdpb_util/Verbosity.hxx"

#include <boost/program_options.hpp>
#include <boost/property_tree/ptree.hpp>

#include <filesystem>

struct Pmp2functions_Parameters
{
  int precision;
  std::filesystem::path input_file;
  std::filesystem::path output_path;
  Verbosity verbosity;

  std::vector<std::string> command_arguments;

  Pmp2functions_Parameters(int argc, char **argv);
  [[nodiscard]] bool is_valid() const;
};

std::ostream &operator<<(std::ostream &os, const Pmp2functions_Parameters &p);

boost::property_tree::ptree
to_property_tree(const Pmp2functions_Parameters &p);