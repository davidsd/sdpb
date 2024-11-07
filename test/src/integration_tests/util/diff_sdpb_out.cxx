#include "diff.hxx"

#include "Float.hxx"
#include "Parse_SDP.hxx"
#include "json.hxx"
#include "sdpb_util/to_matrix.hxx"

#include <catch2/catch_amalgamated.hpp>

#include <rapidjson/document.h>
#include <rapidjson/istreamwrapper.h>

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

  struct Parse_Iterations_Json : boost::noncopyable
  {
    std::vector<std::vector<std::pair<std::string, std::string>>> iterations;
    explicit Parse_Iterations_Json(const fs::path &path)
    {
      CAPTURE(path);
      REQUIRE(exists(path));
      std::ifstream is(path);
      rapidjson::IStreamWrapper wrapper(is);
      rapidjson::Document document;
      document.ParseStream(wrapper);

      REQUIRE(document.IsArray());
      for(const auto &json_iteration : document.GetArray())
        {
          iterations.emplace_back();
          REQUIRE(json_iteration.IsObject());
          const auto obj = json_iteration.GetObject();
          for(auto it = obj.MemberBegin(); it != obj.MemberEnd(); ++it)
            {
              {
                CAPTURE(it->name.GetType());
                REQUIRE(it->name.IsString());
              }
              auto name = std::string(it->name.GetString(),
                                      it->name.GetStringLength());

              std::string value;
              CAPTURE(it->value.GetType());
              if(it->value.IsString())
                {
                  value = std::string(it->value.GetString(),
                                      it->value.GetStringLength());
                }
              else if(it->value.IsInt())
                {
                  value = std::to_string(it->value.GetInt());
                }
              else if(it->value.IsNumber())
                {
                  value = std::to_string(it->value.GetDouble());
                }
              else
                {
                  FAIL("Unsupported JSON value type");
                }
              iterations.back().emplace_back(name, value);
            }
        }
    }
  };

  struct Parse_c_minus_By_Json
  {
    std::vector<Float_Vector> blocks;
    explicit Parse_c_minus_By_Json(const fs::path &path)
    {
      CAPTURE(path);
      REQUIRE(exists(path));
      std::ifstream is(path);
      rapidjson::IStreamWrapper wrapper(is);
      rapidjson::Document document;
      document.ParseStream(wrapper);

      REQUIRE(document.IsObject());
      REQUIRE(document.HasMember("c_minus_By"));
      REQUIRE(document["c_minus_By"].IsArray());
      for(const auto &block : document["c_minus_By"].GetArray())
        {
          blocks.emplace_back(Test_Util::Json::parse_Float_Vector(block));
        }
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

  void diff_iterations_json(const fs::path &a_iterations_json,
                            const fs::path &b_iterations_json)
  {
    CAPTURE(a_iterations_json);
    CAPTURE(b_iterations_json);
    const auto a_iterations
      = Parse_Iterations_Json(a_iterations_json).iterations;
    const auto b_iterations
      = Parse_Iterations_Json(b_iterations_json).iterations;
    DIFF(a_iterations.size(), b_iterations.size());
    for(size_t iter_index = 0; iter_index < a_iterations.size(); ++iter_index)
      {
        const size_t iteration = iter_index + 1;
        CAPTURE(iteration);
        const auto &a_iter = a_iterations.at(iter_index);
        const auto &b_iter = b_iterations.at(iter_index);
        DIFF(a_iter.size(), b_iter.size());
        for(size_t index = 0; index < a_iter.size(); ++index)
          {
            CAPTURE(index);
            const auto [a_key, a_value_str] = a_iter.at(index);
            const auto [b_key, b_value_str] = b_iter.at(index);
            DIFF(a_key, b_key);
            // Timings can vary, we ignore them
            if(a_key == "total_time" || a_key == "iter_time")
              continue;
            // "block_name" can be different if several blocks
            // have the same condition number up to rounding errors.
            // For example, this happened in test/out/dfibo-0-0-j=3-c=3.0000-d=3-s=6/pmp.xml:
            // for iteration=1, block_name="block_2" for format=json and block_name="block_19" for format=bin
            if(a_key == "block_name")
              continue;

            CAPTURE(a_key);
            CAPTURE(a_value_str);
            CAPTURE(b_value_str);

            El::BigFloat a_value = a_value_str;
            El::BigFloat b_value = b_value_str;

            // Errors may differ as they approach zero,
            // so we do not compare small values
            if(a_key == "P-err" || a_key == "p-err" || a_key == "D-err"
               || a_key == "R-err")
              {
                auto prec = Test_Util::REQUIRE_Equal::diff_precision > 0
                              ? Test_Util::REQUIRE_Equal::diff_precision
                              : El::gmp::Precision();
                // Trying to be conservative, otherwise skydiving test fails:
                prec /= 2;
                const auto abs_eps = El::BigFloat(1) >>= prec;
                if(Abs(a_value) + Abs(b_value) < abs_eps)
                  continue;
              }

            DIFF(a_value, b_value);
          }
      }
  }

  void diff_c_minus_By_json(const fs::path &a_json, const fs::path &b_json)
  {
    CAPTURE(a_json);
    CAPTURE(b_json);
    const Parse_c_minus_By_Json a(a_json);
    const Parse_c_minus_By_Json b(b_json);
    DIFF(a.blocks, b.blocks);
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

    const std::string c_minus_By = "c_minus_By";
    const fs::path c_minus_By_rel_path
      = fs::path("c_minus_By") / "c_minus_By.json";

    std::vector<std::string> my_filenames = filenames;
    if(my_filenames.empty())
      {
        fs::directory_iterator a_it(a_out_dir);
        fs::directory_iterator b_it(b_out_dir);
        fs::directory_iterator end{};
        for(const auto &a : fs::directory_iterator(a_out_dir))
          {
            CAPTURE(a);
            // Possible directories:
            // out/c_minus_By/, out/c_minus_By.0/, out/c_minus_By.1/ etc.
            if(a.is_directory())
              {
                const auto dirname = a.path().filename().string();
                REQUIRE(dirname.substr(0, c_minus_By.size()) == c_minus_By);

                // check only the last c_minus_By file, out/c_minus_By/c_minus_By.json
                if(dirname == c_minus_By)
                  {
                    my_filenames.push_back(c_minus_By_rel_path);
                  }
              }
            else
              {
                REQUIRE(is_regular_file(a.path()));
                my_filenames.push_back(a.path().filename().string());
              }
          }
        // Check that all files from b_out_dir exist in a_out_dir
        for(const auto &b : fs::recursive_directory_iterator(b_out_dir))
          {
            CAPTURE(b);
            if(!b.is_regular_file())
              continue;
            const auto rel_path = fs::relative(b.path(), b_out_dir);
            if(rel_path.has_parent_path())
              {
                const auto parent_dir = rel_path.parent_path().string();
                if(parent_dir != c_minus_By)
                  {
                    // Skip out/c_minus_By.0/, out/c_minus_By.1/ etc.
                    REQUIRE(parent_dir.substr(0, c_minus_By.size())
                            == c_minus_By);
                    continue;
                  }
              }
            REQUIRE(is_regular_file(a_out_dir / rel_path));
          }
      }

    for(const auto &name : my_filenames)
      {
        auto a = a_out_dir / name;
        auto b = b_out_dir / name;

        if(name == "out.txt")
          diff_sdpb_out_txt(a, b, out_txt_keys);
        else if(name.substr(0, 10) == "iterations")
          diff_iterations_json(a, b);
        else if(name == c_minus_By_rel_path)
          diff_c_minus_By_json(a, b);
        else
          diff_matrix_txt(a, b);
      }
  }

  void check_c_minus_By(const fs::path &sdp_dir, const fs::path &sdpb_out_dir,
                        const unsigned int input_precision,
                        const unsigned int diff_precision,
                        const Test_Case_Runner &runner)
  {
    Float_Binary_Precision prec(input_precision, diff_precision);

    CAPTURE(sdp_dir);
    CAPTURE(sdpb_out_dir);

    const auto vec_to_matrix
      = [&](const Float_Vector &vec) { return to_matrix(std::vector({vec})); };

    const Parse_SDP sdp(sdp_dir,runner);
    const auto y_vec = Parse_Matrix_Txt(sdpb_out_dir / "y.txt").elements;
    const Float_Matrix y = vec_to_matrix(y_vec);

    const auto c_minus_By
      = Parse_c_minus_By_Json(sdpb_out_dir / "c_minus_By" / "c_minus_By.json")
          .blocks;

    DIFF(sdp.block_data.size(), c_minus_By.size());

    for(size_t block_index = 0; block_index < sdp.block_data.size();
        block_index++)
      {
        CAPTURE(block_index);
        const auto &block_data = sdp.block_data.at(block_index);
        const auto &c = vec_to_matrix(block_data.constraint_constants);
        const auto &B = block_data.constraint_matrix;

        const Float_Matrix c_minus_By_block
          = vec_to_matrix(c_minus_By.at(block_index));

        auto c_minus_By_block_2 = c;
        // c-> - B.y + c
        El::Gemv(El::NORMAL, Float(-1), B, y, Float(1), c_minus_By_block_2);
        DIFF(c_minus_By_block, c_minus_By_block_2);
      }
  }
}
