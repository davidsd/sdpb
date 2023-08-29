#include "integration_tests/common.hxx"

#include <filesystem>

namespace fs = std::filesystem;

TEST_CASE("sdp2input")
{
  INFO("Simple sdp2input tests for different formats: .nsv, .m, .json");

  auto data_dir = Test_Config::test_data_dir / "sdp2input";

  Test_Util::Test_Case_Runner::Named_Args_Map default_args{
    {"--precision", "512"},
    {"--debug", "true"},
  };
  for(std::string format : {"bin", "json"})
    {
      DYNAMIC_SECTION(format)
      {
        auto sdp_orig = data_dir / ("sdp_json.orig.zip");
        for(std::string input_name :
            {"sdp2input_test.json", "sdp2input_split.nsv", "sdp2input_test.m"})
          {
            DYNAMIC_SECTION(input_name)
            {
              Test_Util::Test_Case_Runner runner("sdp2input/" + format + "/"
                                                 + input_name);

              Test_Util::Test_Case_Runner::Named_Args_Map args(default_args);
              args["--input"] = (data_dir / input_name).string();
              auto sdp_zip = (runner.output_dir / "sdp.zip").string();
              args["--output"] = sdp_zip;
              args["--outputFormat"] = format;

              runner.create_nested("run").mpi_run({"build/sdp2input"}, args);

              {
                INFO("Check that sdp2input actually uses --outputFormat="
                     << format);
                auto sdp_unzip
                  = runner.create_nested("format").unzip_to_temp_dir(sdp_zip);
                auto block_data_0_path
                  = sdp_unzip / ("block_data_0." + format);
                CAPTURE(block_data_0_path);
                REQUIRE(is_regular_file(block_data_0_path));
              }

              Test_Util::REQUIRE_Equal::diff_sdp_zip(
                sdp_zip, sdp_orig, 512, 512, runner.create_nested("diff"));

              REQUIRE(fs::file_size(sdp_zip + ".profiling.0") > 0);
              REQUIRE(fs::file_size(sdp_zip + ".profiling.1") > 0);
            }
          }
      }
    }
}
