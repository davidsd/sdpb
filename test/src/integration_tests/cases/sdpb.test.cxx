#include "integration_tests/common.hxx"
#include <boost/filesystem.hpp>

TEST_CASE("sdpb")
{
  auto data_dir = Test_Config::test_data_dir / "sdpb";
  // default sdp input file
  auto sdp_path = Test_Config::test_data_dir / "sdp.zip";

  const Test_Util::Test_Case_Runner::Named_Args_Map default_args{
    {"--sdpDir", sdp_path.string()},
    {"--precision", "1024"},
    {"--noFinalCheckpoint", ""},
    {"--procsPerNode", "1"},
  };

  // Set appropriate checkpoint and out dirs, the run SDPB.
  auto run_sdpb_set_out_ck_dirs
    = [](const Test_Util::Test_Case_Runner &runner,
         Test_Util::Test_Case_Runner::Named_Args_Map &args) -> int {
    args["--checkpointDir"] = (runner.output_dir / "ck").string();
    args["--outDir"] = (runner.output_dir / "out").string();
    return runner.mpi_run({"build/sdpb"}, args);
  };

  // create file with readonly premissions
  auto create_readonly = [](const boost::filesystem::path &path) {
    create_directories(path.parent_path());
    boost::filesystem::ofstream os(path);
    os << "";
    boost::filesystem::permissions(path, boost::filesystem::others_read);
  };

  SECTION("sdpb/sdpb")
  {
    Test_Util::Test_Case_Runner runner("sdpb");
    auto args = default_args;
    int res_run = run_sdpb_set_out_ck_dirs(runner.create_nested("run"), args);
    REQUIRE(res_run == 0);
    int res_diff = runner.create_nested("diff").diff(
      args["--outDir"], data_dir / "test_out_orig");
    REQUIRE(res_diff == 0);
  }

  // Check that sdpb fails on different IO errors
  SECTION("sdpb/sdpb_io_tests")
  {
    {
      Test_Util::Test_Case_Runner runner("sdpb/io_tests/write_profile");
      auto args = default_args;
      args["--maxIterations"] = "1";
      args["--verbosity"] = "2";

      create_readonly(runner.output_dir / "ck.profiling.0");
      int res_run = run_sdpb_set_out_ck_dirs(runner, args);
      REQUIRE(res_run == 1);
      REQUIRE(runner.stderr_contains_substring("Error when writing to"));
      REQUIRE(
        boost::filesystem::file_size(runner.output_dir / "ck.profiling.1")
        > 0);
    }

    for(std::string name :
        {"out.txt", "x_0.txt", "y.txt", "X_matrix_0.txt", "Y_matrix_0.txt"})
      {
        Test_Util::Test_Case_Runner runner("sdpb/io_tests/" + name);
        auto args = default_args;
        args["--maxIterations"] = "1";
        args["--writeSolution"] = "x,y,X,Y";

        create_readonly(runner.output_dir / "out" / name);
        int res_run = run_sdpb_set_out_ck_dirs(runner, args);
        REQUIRE(res_run == 1);
        bool cannot_write
          = runner.stderr_contains_substring("Cannot write to")
            || runner.stderr_contains_substring("Error when writing to");
        REQUIRE(cannot_write);
      }

    {
      Test_Util::Test_Case_Runner runner("sdpb/io_tests/input_corruption");
      auto sdp_corrupted = runner.output_dir / "sdp_corrupted.zip";
      boost::filesystem::create_directories(runner.output_dir);
      boost::filesystem::copy(sdp_path, sdp_corrupted,
                              boost::filesystem::copy_options::recursive);
      boost::filesystem::ofstream os(sdp_corrupted);
      os << "any bytes to corrupt zip archive";

      auto args = default_args;
      args["--sdpDir"] = sdp_corrupted.string();
      args["--maxIterations"] = "1";

      int res_run = run_sdpb_set_out_ck_dirs(runner, args);
      REQUIRE(res_run == 1);
    }

    {
      Test_Util::Test_Case_Runner runner("sdpb/io_tests/checkpoint_read");
      auto args = default_args;
      args["--maxIterations"] = "1";
      args["--writeSolution"] = "x,y,X,Y";

      int res_run = run_sdpb_set_out_ck_dirs(runner, args);
      REQUIRE(res_run == 0);

      // read from checkpoint successfully
      int res_run_ck = run_sdpb_set_out_ck_dirs(runner, args);
      REQUIRE(res_run_ck == 0);

      // now use outDir as checkpoint,
      // remove read permissions => fail to read
      args["--checkpointDir"] = args["--outDir"];
      boost::filesystem::permissions(runner.output_dir / "out"
                                       / "X_matrix_0.txt",
                                     boost::filesystem::perms::no_perms);
      const Test_Util::Test_Case_Runner runner_noread
        = runner.create_nested("noread");
      int res_run_ck_noread = runner_noread.mpi_run({"build/sdpb"}, args);
      REQUIRE(res_run_ck_noread == 1);
      REQUIRE(runner_noread.stderr_contains_substring(
        "Unable to open checkpoint file"));
    }

    {
      Test_Util::Test_Case_Runner runner("sdpb/io_tests/checkpoint_corrupt");
      auto args = default_args;
      args["--maxIterations"] = "1";
      args["--writeSolution"] = "x,y,X,Y";
      int res_run = run_sdpb_set_out_ck_dirs(runner, args);
      REQUIRE(res_run == 0);

      // corrupt X_matrix: remove all data after first two lines
      auto X_matrix
        = (boost::filesystem::path(args["--outDir"]) / "X_matrix_0.txt");
      // remove all but first two lines from X_matrix
      {
        boost::filesystem::ifstream is(X_matrix);
        std::string line_0, line_1;
        std::getline(is, line_0);
        std::getline(is, line_1);
        is.close();

        boost::filesystem::ofstream os(X_matrix);
        os << line_0 << std::endl;
        os << line_1 << std::endl;
      }

      // try to read checkpoint with corrupted X_matrix
      args["--checkpointDir"] = args["--outDir"];
      const Test_Util::Test_Case_Runner runner_corrupt
        = runner.create_nested("read_corrupt");
      int res_run_ck_corrupt = runner_corrupt.mpi_run({"build/sdpb"}, args);
      REQUIRE(res_run_ck_corrupt == 1);
      REQUIRE(
        runner_corrupt.stderr_contains_substring("Corrupted data in file"));
    }
  }
}