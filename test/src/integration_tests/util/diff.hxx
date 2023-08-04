#pragma once

#include "Float.hxx"
#include "Test_Case_Runner.hxx"

#include <catch2/catch_amalgamated.hpp>
#include <boost/filesystem.hpp>

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
  void diff_sdpb_output_dir(const boost::filesystem::path &a_out_dir,
                            const boost::filesystem::path &b_out_dir,
                            unsigned int binary_precision,
                            const std::vector<std::string> &filenames = {},
                            const std::vector<std::string> &out_txt_keys = {});

  void diff_sdp_zip(const boost::filesystem::path &a_sdp_zip,
                    const boost::filesystem::path &b_sdp_zip,
                    unsigned int binary_precision, Test_Case_Runner runner);

  inline void diff(int a, int b) { REQUIRE(a == b); }
  inline void diff(const Float &a, const Float &b)
  {
    // a and b will not be printed with all digits,
    // so we display also their (sometimes small) difference.
    CAPTURE(a - b);
    REQUIRE(a == b);
  }
  template <class T>
  inline void diff(const std::vector<T> &a, const std::vector<T> &b)
  {
    REQUIRE(a.size() == b.size());
    for(size_t i = 0; i < a.size(); ++i)
      {
        CAPTURE(i);
        diff(a[i], b[i]);
      }
  }
}
