#include "integration_tests/common.hxx"

#include <boost/filesystem.hpp>

TEST_CASE("sdp2input")
{
  INFO("Simple sdp2input tests for different formats: .nsv, .m, .json");

  auto data_dir = Test_Config::test_data_dir / "sdp2input";

  Test_Util::Test_Case_Runner::Named_Args_Map default_args{
    {"--precision", "512"},
    {"--debug", "true"},
  };

  auto sdp_orig = data_dir / ("sdp_orig.zip");
  for(std::string input_name :
      {"sdp2input_test.json", "sdp2input_split.nsv", "sdp2input_test.m"})
    {
      DYNAMIC_SECTION(input_name)
      {
        Test_Util::Test_Case_Runner runner("sdp2input/" + input_name);

        Test_Util::Test_Case_Runner::Named_Args_Map args(default_args);
        args["--input"] = (data_dir / input_name).string();
        auto sdp_zip = (runner.output_dir / "sdp.zip").string();
        args["--output"] = sdp_zip;

        runner.create_nested("run").mpi_run({"build/sdp2input"}, args);

        Test_Util::REQUIRE_Equal::diff_sdp_zip(sdp_zip, sdp_orig, 1024,
                                               runner.create_nested("diff"));

        REQUIRE(boost::filesystem::file_size(sdp_zip + ".profiling.0") > 0);
        REQUIRE(boost::filesystem::file_size(sdp_zip + ".profiling.1") > 0);
      }
    }
}
