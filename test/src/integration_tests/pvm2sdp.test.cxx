#include "common.hxx"

#include <boost/filesystem.hpp>

TEST_CASE("pvm2sdp")
{
  Test_Util::Test_Case_Runner main_runner("pvm2sdp");
  auto data_dir = main_runner.data_dir;
  auto output_dir = main_runner.output_dir;

  auto input = (data_dir / "file_list.nsv").string();
  auto output = (output_dir / "sdp.zip").string();
  int res_run = main_runner.create_nested("run").mpi_run(
    {"build/pvm2sdp 1024", input, output});
  REQUIRE(res_run == 0);

  auto output_orig = (Test_Config::test_data_dir / "sdp.zip").string();
  int res_diff
    = main_runner.create_nested("diff").diff_sdp_zip(output, output_orig);
  REQUIRE(res_diff == 0);

  // prohibit output sdp.zip from writing and check that pvm2sdp fails
  {
    const Test_Util::Test_Case_Runner runner
      = main_runner.create_nested("cannot_write_zip");
    boost::filesystem::create_directories(runner.output_dir);

    auto sdp_readonly_zip = runner.output_dir / "sdp.readonly.zip";
    boost::filesystem::ofstream os(sdp_readonly_zip);
    os << "INVALID ZIP";
    permissions(sdp_readonly_zip, boost::filesystem::others_read);

    int res_nowrite = runner.mpi_run(
      {"build/pvm2sdp 1024", input, sdp_readonly_zip.string()});
    REQUIRE(res_nowrite == 1);
    REQUIRE(runner.stderr_contains_substring(
      "Unable to set options for writing an archive"));
  }

  // pvm2sdp should fail on an incorrect .nsv
  {
    const Test_Util::Test_Case_Runner runner
      = main_runner.create_nested("invalid_nsv");
    auto invalid_nsv = (output_dir / "invalid_file_list.nsv").string();
    {
      boost::filesystem::ofstream os(invalid_nsv);
      os << "no_such_file.xml";
    }

    int res_nsv = runner.mpi_run(
      {"build/pvm2sdp 1024", invalid_nsv, output + ".invalid"});
    REQUIRE(res_nsv == 1);
    REQUIRE(runner.stderr_contains_substring("Unable to parse input file"));
  }
}
