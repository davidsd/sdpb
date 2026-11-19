#pragma once

#include "../Block_File_Format.hxx"
#include "sdpb_util/Verbosity.hxx"

#include <boost/program_options.hpp>
#include <boost/property_tree/ptree.hpp>

#include <filesystem>

struct Pmp2sdp_Parameters
{
  int precision;
  int64_t max_num_poles;
  // Two-set encoding (--encodingPoints): "prefactorMatched" (default: the
  // SDP is written at the approximation points, as without the option),
  // "functionalMatched", "chebyshev" or "file"
  std::string encoding_points = "prefactorMatched";
  double encoding_points_scale = 1;
  std::filesystem::path encoding_scale_file;
  std::filesystem::path encoding_points_file;
  int64_t encoding_extra_precision = 384;
  [[nodiscard]] bool two_set_encoding() const
  {
    return encoding_points != "prefactorMatched";
  }
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