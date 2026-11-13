#include "integration_tests/common.hxx"
#include "integration_tests/util/Parse_SDP.hxx"

#include <filesystem>
#include <map>

namespace fs = std::filesystem;
using namespace std::string_literals;

TEST_CASE("pmp2sdp")
{
  unsigned int precision = 512;
  unsigned int diff_precision = 392;
  Test_Util::Test_Case_Runner::Named_Args_Map default_args{
    {"--precision", std::to_string(precision)},
    {"--verbosity", "debug"},
  };

  SECTION("run")
  {
    std::string input_format = GENERATE("json", "m");
    DYNAMIC_SECTION(input_format)
    {
      auto data_dir = Test_Config::test_data_dir / "pmp2sdp" / input_format;
      auto sdp_orig = data_dir / "sdp_orig";

      for(const std::string &input_filename :
          {"pmp."s + input_format, "file_list.nsv"s})
        DYNAMIC_SECTION(input_filename)
        {
          Test_Util::Test_Case_Runner runner("pmp2sdp/"s + input_format + "/"
                                             + input_filename);
          Test_Util::Test_Case_Runner::Named_Args_Map args(default_args);
          args["--input"] = (data_dir / input_filename).string();
          auto sdp_path = (runner.output_dir / "sdp").string();
          args["--output"] = sdp_path;

          runner.create_nested("run").mpi_run(
            {Test_Config::program("pmp2sdp"), args});

          Test_Util::REQUIRE_Equal::diff_sdp(sdp_path, sdp_orig, precision,
                                             diff_precision,
                                             runner.create_nested("diff"));

          REQUIRE(fs::file_size(sdp_path + ".profiling/profiling.0") > 0);
          REQUIRE(fs::file_size(sdp_path + ".profiling/profiling.1") > 0);
        }
    }
  }

  SECTION("outputFormat")
  {
    INFO("Check different --outputFormat options");
    // TODO use input from "1d" and update output
    auto data_dir
      = Test_Config::test_data_dir / "end-to-end_tests" / "1d-old-sampling";

    for(std::string output_format : {"", "bin", "json"})
      {
        auto format_description
          = output_format.empty() ? "default(bin)" : output_format;
        CAPTURE(output_format);
        CAPTURE(format_description);

        DYNAMIC_SECTION(format_description)
        {
          Test_Util::Test_Case_Runner runner("pmp2sdp/outputFormat"s + "/"
                                             + format_description);
          Test_Util::Test_Case_Runner::Named_Args_Map args(default_args);
          args["--input"] = (data_dir / "input" / "pmp.json").string();
          auto sdp_path = runner.output_dir / "sdp";
          args["--output"] = sdp_path.string();
          if(!output_format.empty())
            args["--outputFormat"] = output_format;
          runner.create_nested("run").mpi_run(
            {Test_Config::program("pmp2sdp"), args});

          {
            INFO("Check that pmp2sdp actually uses --outputFormat="
                 << format_description);
            auto block_data_0_path
              = sdp_path
                / ("block_data_0."
                   + (output_format.empty() ? "bin" : output_format));
            CAPTURE(block_data_0_path);
            REQUIRE(is_regular_file(block_data_0_path));
          }

          auto sdp_orig = data_dir / "output" / "sdp";

          Test_Util::REQUIRE_Equal::diff_sdp(sdp_path, sdp_orig, precision,
                                             diff_precision,
                                             runner.create_nested("diff"));
        }
      }
  }

  SECTION("encodingPoints")
  {
    INFO("--encodingPoints (two-set encoding): blocks with more than one "
         "sample point are written at the encoding points; degree-0 blocks "
         "(one sample point) keep the default encoding, and the objective "
         "is unchanged.");
    // dfibo: 98 blocks, 48 of degree 0 (with poles) and 50 with four sample
    // points.  two_set/pmp.json: a scalar block of degree 6 and a 2x2 block
    // of degree 4 (with poles), a constant block without poles (its one
    // sample point is x = 0) and a degree-1 block.
    struct Input
    {
      std::string name;
      fs::path path;
      size_t num_degree_0, num_reencoded;
    };
    const Input input = GENERATE(
      Input{"dfibo",
            Test_Config::test_data_dir / "end-to-end_tests"
              / "dfibo-0-0-j=3-c=3.0000-d=3-s=6" / "input" / "pmp.xml",
            48, 50},
      Input{"two_set", Test_Config::test_data_dir / "pmp2sdp" / "two_set"
                         / "pmp.json",
            1, 3});
    CAPTURE(input.name);
    Test_Util::Test_Case_Runner runner("pmp2sdp/encodingPoints/"
                                       + input.name);
    std::map<std::string, fs::path> sdp_paths;
    for(const std::string rule :
        {"prefactorMatched", "functionalMatched", "chebyshev"})
      {
        Test_Util::Test_Case_Runner::Named_Args_Map args(default_args);
        args["--input"] = input.path.string();
        sdp_paths[rule] = runner.output_dir / (rule + ".sdp");
        args["--output"] = sdp_paths[rule].string();
        args["--outputFormat"] = "json";
        args["--encodingPoints"] = rule;
        runner.create_nested(rule).mpi_run(
          {Test_Config::program("pmp2sdp"), args});
      }

    using namespace Test_Util::REQUIRE_Equal;
    Float_Binary_Precision prec(precision, diff_precision);
    const Parse_SDP sdp_default(sdp_paths.at("prefactorMatched"),
                                runner.create_nested("parse-default"));
    for(const std::string rule : {"functionalMatched", "chebyshev"})
      DYNAMIC_SECTION(rule)
      {
        const Parse_SDP sdp(sdp_paths.at(rule),
                            runner.create_nested("parse-" + rule));
        DIFF(sdp.control.num_blocks, sdp_default.control.num_blocks);
        DIFF(sdp.objectives.constant, sdp_default.objectives.constant);
        DIFF(sdp.objectives.b, sdp_default.objectives.b);
        size_t num_degree_0 = 0, num_reencoded = 0;
        for(int block = 0; block < sdp.control.num_blocks; ++block)
          {
            CAPTURE(block);
            const auto &info = sdp.block_info.at(block);
            const auto &data = sdp.block_data.at(block);
            const auto &data_default = sdp_default.block_data.at(block);
            DIFF(info.dim, sdp_default.block_info.at(block).dim);
            DIFF(info.num_points, sdp_default.block_info.at(block).num_points);
            if(info.num_points == 1)
              {
                ++num_degree_0;
                DIFF(data.bilinear_bases_even,
                     data_default.bilinear_bases_even);
                DIFF(data.bilinear_bases_odd, data_default.bilinear_bases_odd);
                DIFF(data.constraint_constants,
                     data_default.constraint_constants);
                DIFF(data.constraint_matrix, data_default.constraint_matrix);
              }
            else
              {
                // the rows are at other points: c differs from the default's
                ++num_reencoded;
                REQUIRE(data.constraint_constants.size()
                        == data_default.constraint_constants.size());
                bool differs = false;
                for(size_t i = 0; i < data.constraint_constants.size(); ++i)
                  {
                    const auto &x = data.constraint_constants.at(i);
                    const auto &y = data_default.constraint_constants.at(i);
                    if(El::Abs(x - y) > El::BigFloat("1e-20")
                                          * (El::Abs(x) + El::Abs(y)))
                      differs = true;
                  }
                REQUIRE(differs);
              }
          }
        REQUIRE(num_degree_0 == input.num_degree_0);
        REQUIRE(num_reencoded == input.num_reencoded);
      }
  }

  SECTION("filesystem errors")
  {
    INFO("pmp2sdp should fail due to invalid input/output arguments");
    auto data_dir = Test_Config::test_data_dir / "end-to-end_tests" / "1d";
    int num_procs = 2;

    auto input = (data_dir / "input" / "pmp.json").string();

    SECTION("sdp_temp_not_empty")
    {
      const Test_Util::Test_Case_Runner runner("pmp2sdp/sdp_temp_not_empty");
      fs::create_directories(runner.output_dir);

      auto sdp_dir = runner.output_dir / "sdp";

      {
        auto sdp_temp_dir = runner.output_dir / "sdp_temp";
        fs::create_directories(runner.output_dir / "sdp_temp");
        std::ofstream os(sdp_temp_dir / "some_file.json");
        os << "some data";
      }

      Test_Util::Test_Case_Runner::Named_Args_Map args(default_args);
      args["--input"] = input;
      args["--output"] = sdp_dir.string();
      // We allow pmp2sdp to overwrite exsiting sdp_temp
      runner.mpi_run({Test_Config::program("pmp2sdp"), args}, num_procs, 0,
                     "exists and will be overwritten");
    }
    SECTION("sdp_dir_exists")
    {
      INFO("Prohibit writing to existing directory:");
      const Test_Util::Test_Case_Runner runner("pmp2sdp/sdp_dir_exists");
      fs::create_directories(runner.output_dir);

      auto sdp_dir = runner.output_dir / "sdp";
      {
        fs::create_directories(sdp_dir);
        std::ofstream os(sdp_dir / "control.json");
        os << "some data";
      }

      Test_Util::Test_Case_Runner::Named_Args_Map args(default_args);
      args["--input"] = input;
      args["--output"] = sdp_dir.string();
      // We allow pmp2sdp to overwrite exsiting sdp
      runner.mpi_run({Test_Config::program("pmp2sdp"), args}, num_procs, 0,
                     "exists and will be overwritten");
    }

    SECTION("cannot_write_zip")
    {
      INFO("Prohibit output sdp.zip from writing and check that pmp2sdp "
           "fails:");
      const Test_Util::Test_Case_Runner runner("pmp2sdp/cannot_write_zip");
      fs::create_directories(runner.output_dir);

      auto sdp_readonly_zip = runner.output_dir / "sdp.readonly.zip";
      std::ofstream os(sdp_readonly_zip);
      os << "INVALID ZIP";
      permissions(sdp_readonly_zip, fs::perms::others_read);

      Test_Util::Test_Case_Runner::Named_Args_Map args(default_args);
      args["--input"] = input;
      args["--output"] = sdp_readonly_zip.string();
      args["--zip"] = "";
      // We allow pmp2sdp to overwrite exsiting sdp
      runner.mpi_run({Test_Config::program("pmp2sdp"), args}, num_procs, 0,
                     "exists and will be overwritten");
    }

    SECTION("invalid_nsv")
    {
      INFO("pmp2sdp should fail on an incorrect .nsv");
      const Test_Util::Test_Case_Runner runner("pmp2sdp/invalid_nsv");
      fs::create_directories(runner.output_dir);
      auto invalid_nsv
        = (runner.output_dir / "invalid_file_list.nsv").string();
      {
        std::ofstream os(invalid_nsv);
        os << "no_such_file.xml";
      }

      Test_Util::Test_Case_Runner::Named_Args_Map args(default_args);
      args["--input"] = invalid_nsv;
      auto sdp_path = (runner.output_dir / "sdp").string();
      args["--output"] = sdp_path;

      auto sdp_invalid_zip = (runner.output_dir / "sdp.invalid.zip").string();
      runner.mpi_run({Test_Config::program("pmp2sdp"), args}, num_procs, 1,
                     "No such file");
    }
  }

  SECTION("duplicate objectives")
  {
    INFO("Two PMP files having same objective vector.");
    INFO("pmp2sdp should not fail, see "
         "https://github.com/davidsd/sdpb/issues/251");
    INFO("Same input as end-to-end_tests/1d-constraints/pmp.xml, but split "
         "into two files");

    auto data_dir
      = Test_Config::test_data_dir / "pmp2sdp" / "duplicate-objectives";
    int num_procs = GENERATE(1, 2);
    auto input_nsv = (data_dir / "pmp.nsv").string();

    DYNAMIC_SECTION("num_procs=" << num_procs)
    {
      const Test_Util::Test_Case_Runner runner(
        "pmp2sdp/duplicate-objectives/num_procs=" + std::to_string(num_procs));

      Test_Util::Test_Case_Runner::Named_Args_Map args(default_args);
      args["--input"] = input_nsv;
      auto sdp_path = (runner.output_dir / "sdp").string();
      args["--output"] = sdp_path;

      runner.create_nested("run").mpi_run(
        {Test_Config::program("pmp2sdp"), args}, num_procs);

      auto sdp_orig = Test_Config::test_data_dir / "end-to-end_tests"
                      / "1d-constraints" / "output" / "sdp";
      Test_Util::REQUIRE_Equal::diff_sdp(sdp_path, sdp_orig, precision,
                                         diff_precision,
                                         runner.create_nested("diff"));
    }
  }

  SECTION("objectives are different")
  {
    INFO("pmp2sdp should fail because two files have different objectives");

    auto data_dir
      = Test_Config::test_data_dir / "pmp2sdp" / "wrong-objectives";
    int num_procs = GENERATE(1, 2);
    auto input_nsv = (data_dir / "pmp.nsv").string();

    DYNAMIC_SECTION("num_procs=" << num_procs)
    {
      const Test_Util::Test_Case_Runner runner(
        "pmp2sdp/wrong-objectives/num_procs=" + std::to_string(num_procs));

      Test_Util::Test_Case_Runner::Named_Args_Map args(default_args);
      args["--input"] = input_nsv;
      auto sdp_path = (runner.output_dir / "sdp").string();
      args["--output"] = sdp_path;
      runner.mpi_run({Test_Config::program("pmp2sdp"), args}, num_procs, 1,
                     "Found different objective vectors in input files:");
    }
  }
}
