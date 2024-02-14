#include "diff.hxx"

#include "Float.hxx"

#include <catch2/catch_amalgamated.hpp>

#include <boost/algorithm/algorithm.hpp>
#include <boost/algorithm/string.hpp>

namespace fs = std::filesystem;

// Parser classes
namespace
{
  struct Parse_Matrix_Txt
  {
    int height;
    int width;
    Float_Vector elements;
    explicit Parse_Matrix_Txt(const fs::path &path)
    {
      CAPTURE(path);
      REQUIRE(is_regular_file(path));
      std::ifstream is(path);
      height = width = 0;
      is >> height;
      is >> width;
      REQUIRE(height > 0);
      REQUIRE(width > 0);
      elements.reserve(height * width);
      for(int i = 0; i < height; ++i)
        {
          for(int k = 0; k < width; ++k)
            {
              Float f;
              is >> f;
              elements.push_back(f);
            }
        }
    }
  };

  // out.txt produced by SDPB
  struct Parse_Sdpb_Out_Txt : boost::noncopyable
  {
    std::string terminate_reason;
    std::map<std::string, Float> float_map;

    explicit Parse_Sdpb_Out_Txt(const fs::path &path)
    {
      CAPTURE(path);
      REQUIRE(is_regular_file(path));
      std::ifstream is(path);
      std::string line;
      while(std::getline(is, line))
        {
          // line format:
          // key = value;
          CAPTURE(line);
          std::vector<std::string> tokens;
          boost::split(tokens, line, boost::is_any_of("=;"));
          if(tokens.size() <= 1)
            break; // empty string (the last one)

          REQUIRE(tokens.size() == 3); // key, value, empty "" at the end
          auto key = tokens[0];
          auto value = tokens[1];
          boost::trim(key);
          boost::trim(value);
          if(key == "terminateReason")
            {
              terminate_reason = value;
            }
          else
            {
              float_map[key] = Float(value);
            }
        }
      REQUIRE(!terminate_reason.empty());
    }
  };
}

// Helper functions
namespace
{
  using Test_Util::REQUIRE_Equal::diff;
  // Compare matrices written by SDPB save_solution() method,
  // e.g. y.txt or X.txt
  void
  diff_matrix_txt(const fs::path &a_matrix_txt, const fs::path &b_matrix_txt)
  {
    CAPTURE(a_matrix_txt);
    CAPTURE(b_matrix_txt);
    Parse_Matrix_Txt a(a_matrix_txt);
    Parse_Matrix_Txt b(b_matrix_txt);
    DIFF(a.height, b.height);
    DIFF(a.width, b.width);
    DIFF(a.elements, b.elements);
  }

  // Compare out.txt
  void diff_sdpb_out_txt(const fs::path &a_out_txt, const fs::path &b_out_txt,
                         const std::vector<std::string> &keys_to_compare)
  {
    CAPTURE(a_out_txt);
    CAPTURE(b_out_txt);
    Parse_Sdpb_Out_Txt a(a_out_txt);
    Parse_Sdpb_Out_Txt b(a_out_txt);

    auto keys = keys_to_compare;
    // By default, test each key except for "Solver runtime"
    if(keys.empty())
      keys = {"terminateReason", "primalObjective", "dualObjective",
              "dualityGap",      "primalError",     "dualError"};
    CAPTURE(keys);

    for(const auto &key : keys)
      {
        CAPTURE(key);
        if(key == "terminateReason")
          {
            REQUIRE(a.terminate_reason == b.terminate_reason);
            continue;
          }
        if(key == "Solver runtime")
          WARN("Solver runtime may differ for different runs, "
               "do you really want to check it?");

        auto &a_map = a.float_map;
        auto &b_map = b.float_map;
        REQUIRE(a_map.find(key) != a_map.end());
        REQUIRE(b_map.find(key) != b_map.end());

        DIFF(a_map[key], b_map[key]);
      }
  }
}

// Implementation
namespace Test_Util::REQUIRE_Equal
{
  void
  diff_sdpb_output_dir(const fs::path &a_out_dir, const fs::path &b_out_dir,
                       unsigned int input_precision,
                       unsigned int diff_precision,
                       const std::vector<std::string> &filenames,
                       const std::vector<std::string> &out_txt_keys)
  {
    INFO("diff sdpb output");
    REQUIRE(a_out_dir != b_out_dir);
    CAPTURE(a_out_dir);
    CAPTURE(b_out_dir);
    REQUIRE(is_directory(a_out_dir));
    REQUIRE(is_directory(b_out_dir));
    Float_Binary_Precision prec(input_precision, diff_precision);

    std::vector<std::string> my_filenames = filenames;
    if(my_filenames.empty())
      {
        fs::directory_iterator a_it(a_out_dir);
        fs::directory_iterator b_it(b_out_dir);
        fs::directory_iterator end{};
        for(const auto &a : fs::directory_iterator(a_out_dir))
          {
            CAPTURE(a);
            REQUIRE(is_regular_file(a.path()));
            my_filenames.push_back(a.path().filename().string());
          }
        // Check that all files from b_out_dir exist in a_out_dir
        for(const auto &b : fs::directory_iterator(b_out_dir))
          {
            CAPTURE(b);
            REQUIRE(is_regular_file(a_out_dir / b.path().filename()));
          }
      }

    for(const auto &name : my_filenames)
      {
        auto a = a_out_dir / name;
        auto b = b_out_dir / name;

        if(name == "out.txt")
          diff_sdpb_out_txt(a, b, out_txt_keys);
        else
          diff_matrix_txt(a, b);
      }
  }
}
