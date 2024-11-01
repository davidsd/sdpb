#include "integration_tests/common.hxx"

#include <filesystem>

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

          runner.create_nested("run").mpi_run({"build/pmp2sdp"}, args);

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
          runner.create_nested("run").mpi_run({"build/pmp2sdp"}, args);

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
      runner.mpi_run({"build/pmp2sdp"}, args, num_procs, 0,
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
      runner.mpi_run({"build/pmp2sdp"}, args, num_procs, 0,
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
      runner.mpi_run({"build/pmp2sdp"}, args, num_procs, 0,
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
      runner.mpi_run({"build/pmp2sdp"}, args, num_procs, 1, "No such file");
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

      runner.create_nested("run").mpi_run({"build/pmp2sdp"}, args, num_procs);

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
      runner.mpi_run({"build/pmp2sdp"}, args, num_procs, 1,
                     "Found different objective vectors in input files:");
    }
  }
}
