#include "integration_tests/common.hxx"

#include <filesystem>

namespace fs = std::filesystem;

TEST_CASE("pvm2sdp")
{
  INFO("Simple pvm2sdp tests for a one-dimensional problem, same as in "
       "TEST_CASE(sdpb)");
  auto data_dir = Test_Config::test_data_dir / "pvm2sdp";
  int num_procs = 2;

  auto input = (data_dir / "file_list.nsv").string();

  SECTION("pvm2sdp")
  {
    INFO("run pvm2sdp and check output");
    Test_Util::Test_Case_Runner runner("pvm2sdp");
    auto output = (runner.output_dir / "sdp.zip").string();
    auto output_orig = (Test_Config::test_data_dir / "sdp.zip").string();

    for(std::string sdp_format : {"", "bin", "json"})
      {
        // for sdp_format="" pvm2sdp will use bin format (by default)
        auto format_description = sdp_format.empty() ? "default(bin)" : sdp_format;
        DYNAMIC_SECTION(format_description)
        {
          runner.create_nested(format_description + ".run")
            .mpi_run({"build/pvm2sdp", sdp_format, "1024", input, output}, {},
                     2);

          {
            INFO("Check that pvm2sdp actually uses correct outputFormat="
                 << format_description);
            auto sdp_unzip
              = runner.create_nested("format").unzip_to_temp_dir(output);
            auto format = sdp_format.empty() ? "bin" : sdp_format;
            auto block_data_0_path = sdp_unzip / ("block_data_0." + format);
            CAPTURE(block_data_0_path);
            REQUIRE(is_regular_file(block_data_0_path));
          }

          Test_Util::REQUIRE_Equal::diff_sdp_zip(
            output, output_orig, 1024, 1024, runner.create_nested("diff"));
        }
      }
  }

  SECTION("cannot_write_zip")
  {
    INFO("Prohibit output sdp.zip from writing and check that pvm2sdp fails:");
    const Test_Util::Test_Case_Runner runner("pvm2sdp/cannot_write_zip");
    fs::create_directories(runner.output_dir);

    auto sdp_readonly_zip = runner.output_dir / "sdp.readonly.zip";
    std::ofstream os(sdp_readonly_zip);
    os << "INVALID ZIP";
    permissions(sdp_readonly_zip, fs::perms::others_read);

    runner.mpi_run({"build/pvm2sdp 1024", input, sdp_readonly_zip.string()},
                   {}, num_procs, 1,
                   "Unable to set options for writing an archive");
  }

  SECTION("invalid_nsv")
  {
    INFO("pvm2sdp should fail on an incorrect .nsv");
    const Test_Util::Test_Case_Runner runner("pvm2sdp/invalid_nsv");
    fs::create_directories(runner.output_dir);
    auto invalid_nsv = (runner.output_dir / "invalid_file_list.nsv").string();
    {
      std::ofstream os(invalid_nsv);
      os << "no_such_file.xml";
    }

    auto sdp_invalid_zip = (runner.output_dir / "sdp.invalid.zip").string();
    runner.mpi_run({"build/pvm2sdp 1024", invalid_nsv, sdp_invalid_zip}, {},
                   num_procs, 1, "No such file");
  }
}
