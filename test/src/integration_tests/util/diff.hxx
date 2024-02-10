#pragma once

#include "test_util/diff.hxx"
#include "Test_Case_Runner.hxx"

#include <catch2/catch_amalgamated.hpp>
#include <filesystem>

// Functions that REQUIRE equality of given files or folders.
// NB: we call REQUIRE() inside these functions,
// because otherwise useful CAPTURE() information will be lost
namespace Test_Util::REQUIRE_Equal
{
  // Compare SDPB output directories by content.
  // filenames: which files to compare (by default - all)
  // out_txt_keys: which keys compare in out.txt
  // (by default - all except for "Solver runtime")
  // Floating-point numbers are rounded to binary_precision
  void diff_sdpb_output_dir(const std::filesystem::path &a_out_dir,
                            const std::filesystem::path &b_out_dir,
                            unsigned int input_precision,
                            unsigned int diff_precision,
                            const std::vector<std::string> &filenames = {},
                            const std::vector<std::string> &out_txt_keys = {});

  void diff_sdp(const std::filesystem::path &a_sdp,
                const std::filesystem::path &b_sdp,
                unsigned int input_precision, unsigned int diff_precision,
                Test_Case_Runner runner, bool check_normalization = true);

  void diff_functions_json(const std::filesystem::path &a_json,
                           const std::filesystem::path &b_json,
                           unsigned int input_precision,
                           unsigned int diff_precision);

  void
  diff_outer_limits(const std::filesystem::path &a_json,
                    const std::filesystem::path &b_json,
                    unsigned int input_precision, unsigned int diff_precision);

  void
  diff_spectrum(const std::filesystem::path &a_json,
                const std::filesystem::path &b_json,
                unsigned int input_precision, unsigned int diff_precision);
}
