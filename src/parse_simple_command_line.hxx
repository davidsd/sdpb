#pragma once

#include <boost/filesystem.hpp>
#include <vector>

void parse_simple_command_line(
  const std::string &program_name, int argc, char **argv, int &precision,
  std::vector<boost::filesystem::path> &input_files,
  boost::filesystem::path &output_dir);
