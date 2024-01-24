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
    std::string input_format = GENERATE("json", "m", "xml");
    DYNAMIC_SECTION(input_format)
    {
      auto data_dir = Test_Config::test_data_dir / "pmp2sdp" / input_format;
      auto sdp_orig = data_dir / "sdp_orig.zip";

      for(const std::string &input_filename :
          {"pmp."s + input_format, "file_list.nsv"s})
        DYNAMIC_SECTION(input_filename)
        {
          Test_Util::Test_Case_Runner runner("pmp2sdp/"s + input_format + "/"
                                             + input_filename);
          Test_Util::Test_Case_Runner::Named_Args_Map args(default_args);
          args["--input"] = (data_dir / input_filename).string();
          auto sdp_zip = (runner.output_dir / "sdp.zip").string();
          args["--output"] = sdp_zip;

          runner.create_nested("run").mpi_run({"build/pmp2sdp"}, args);

          Test_Util::REQUIRE_Equal::diff_sdp_zip(sdp_zip, sdp_orig, precision,
                                                 diff_precision,
                                                 runner.create_nested("diff"));

          REQUIRE(fs::file_size(sdp_zip + ".profiling/profiling.0") > 0);
          REQUIRE(fs::file_size(sdp_zip + ".profiling/profiling.1") > 0);
        }
    }
  }

  SECTION("outputFormat")
  {
    INFO("Check different --outputFormat options");
    auto data_dir = Test_Config::test_data_dir / "pmp2sdp" / "xml";

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
          args["--input"] = (data_dir / "file_list.nsv").string();
          auto sdp_zip = (runner.output_dir / "sdp.zip").string();
          args["--output"] = sdp_zip;
          if(!output_format.empty())
            args["--outputFormat"] = output_format;
          runner.create_nested("run").mpi_run({"build/pmp2sdp"}, args);

          {
            INFO("Check that pmp2sdp actually uses --outputFormat="
                 << format_description);
            auto sdp_unzip
              = runner.create_nested("unzip").unzip_to_temp_dir(sdp_zip);
            auto block_data_0_path
              = sdp_unzip
                / ("block_data_0."
                   + (output_format.empty() ? "bin" : output_format));
            CAPTURE(block_data_0_path);
            REQUIRE(is_regular_file(block_data_0_path));
          }

          auto sdp_orig = data_dir / "sdp_orig.zip";

          Test_Util::REQUIRE_Equal::diff_sdp_zip(sdp_zip, sdp_orig, precision,
                                                 diff_precision,
                                                 runner.create_nested("diff"));
        }
      }
  }

  SECTION("filesystem errors")
  {
    INFO("pmp2sdp should fail due to invalid input/output arguments");
    auto data_dir = Test_Config::test_data_dir / "pmp2sdp" / "xml";
    int num_procs = 2;

    auto input = (data_dir / "file_list.nsv").string();

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
      auto sdp_zip = (runner.output_dir / "sdp.zip").string();
      args["--output"] = sdp_readonly_zip.string();
      runner.mpi_run({"build/pmp2sdp"}, args, num_procs, 1,
                     "Unable to set options for writing an archive");
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
      auto sdp_zip = (runner.output_dir / "sdp.zip").string();
      args["--output"] = sdp_zip;

      auto sdp_invalid_zip = (runner.output_dir / "sdp.invalid.zip").string();
      runner.mpi_run({"build/pmp2sdp"}, args, num_procs, 1, "No such file");
    }
  }
}
