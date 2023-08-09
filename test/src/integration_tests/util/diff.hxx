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
  inline unsigned int diff_precision;

  // Compare SDPB output directories by content.
  // filenames: which files to compare (by default - all)
  // out_txt_keys: which keys compare in out.txt
  // (by default - all except for "Solver runtime")
  // Floating-point numbers are rounded to binary_precision
  void diff_sdpb_output_dir(const boost::filesystem::path &a_out_dir,
                            const boost::filesystem::path &b_out_dir,
                            unsigned int input_precision,
                            unsigned int diff_precision,
                            const std::vector<std::string> &filenames = {},
                            const std::vector<std::string> &out_txt_keys = {});

  void diff_sdp_zip(const boost::filesystem::path &a_sdp_zip,
                    const boost::filesystem::path &b_sdp_zip,
                    unsigned int input_precision, unsigned int diff_precision,
                    Test_Case_Runner runner);

  void
  diff_outer_limits(const boost::filesystem::path &a_json,
                    const boost::filesystem::path &b_json,
                    unsigned int input_precision, unsigned int diff_precision);

  void
  diff_spectrum(const boost::filesystem::path &a_json,
                const boost::filesystem::path &b_json,
                unsigned int input_precision, unsigned int diff_precision);

  inline void diff(int a, int b) { REQUIRE(a == b); }
  inline void diff(const Float &a, const Float &b)
  {
    if(a == b)
      return;

    CAPTURE(a);
    CAPTURE(b);
    CAPTURE(a - b);

    CAPTURE(diff_precision);
    REQUIRE(diff_precision > 0);

    auto eps = Float(1) >>= diff_precision; // 2^{-precision}
    CAPTURE(eps);
    REQUIRE(Abs(a - b) < eps * (Abs(a) + Abs(b)));
  }
  inline void diff(const std::string &a, const std::string &b)
  {
    REQUIRE(a == b);
  }
  template <class T1, class T2>
  inline void diff(const std::pair<T1, T2> &a, const std::pair<T1, T2> &b)
  {
    INFO("diff std::pair");
    diff(a.first, b.first);
    diff(a.second, b.second);
  }
  template <class T>
  inline void diff(const std::vector<T> &a, const std::vector<T> &b)
  {
    INFO("diff std::vector");
    REQUIRE(a.size() == b.size());
    for(size_t i = 0; i < a.size(); ++i)
      {
        CAPTURE(i);
        diff(a[i], b[i]);
      }
  }
  template <class T>
  inline void diff(const El::Matrix<T> &a, const El::Matrix<T> &b)
  {
    INFO("diff El::Matrix");
    REQUIRE(a.Height() == b.Height());
    REQUIRE(a.Width() == b.Width());
    for(int row = 0; row < a.Height(); ++row)
      for(int col = 0; col < a.Width(); ++col)
        {
          CAPTURE(row);
          CAPTURE(col);
          diff(a.Get(row, col), b.Get(row, col));
        }
  }
}
