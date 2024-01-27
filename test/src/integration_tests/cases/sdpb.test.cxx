#include "integration_tests/common.hxx"
#include <filesystem>

namespace fs = std::filesystem;

TEST_CASE("sdpb")
{
  INFO("SDPB tests for a simple one-dimensional problem:");
  INFO("maximize (-y) s.t. (1 + x^4 + y * (x^4 / 12 + x^2)) >= 0)");
  auto data_dir = Test_Config::test_data_dir / "sdpb";
  // default sdp input file
  auto sdp_path
    = Test_Config::test_data_dir / "pmp2sdp" / "xml" / "sdp_orig.zip";
  CAPTURE(sdp_path);

  int num_procs = 2;

  const Test_Util::Test_Case_Runner::Named_Args_Map default_args{
    {"--sdpDir", sdp_path.string()},
    {"--precision", "1024"},
    {"--noFinalCheckpoint", ""},
  };

  // Set appropriate checkpoint and out dirs, the run SDPB.
  auto run_sdpb_set_out_ck_dirs
    = [](const Test_Util::Test_Case_Runner &runner,
         Test_Util::Test_Case_Runner::Named_Args_Map &args, int num_procs = 2,
         int required_exit_code = 0,
         const std::string &required_error_msg = "") -> void {
    args["--checkpointDir"] = (runner.output_dir / "ck").string();
    args["--outDir"] = (runner.output_dir / "out").string();
    // We removed obsolete --procsPerNode from end-to-end.test.cxx,
    // but we keep it here to check backward compatibility,
    // i.e. SDPB shouldn't fail when we pass this option.
    args["--procsPerNode"] = std::to_string(num_procs);
    runner.mpi_run({"build/sdpb"}, args, num_procs, required_exit_code,
                   required_error_msg);
  };

  // create file with readonly premissions
  auto create_readonly = [](const fs::path &path) {
    create_directories(path.parent_path());
    std::ofstream os(path);
    os << "";
    fs::permissions(path, fs::perms::others_read);
  };

  SECTION("sdpb")
  {
    INFO("Simple sdpb run");
    int sdpb_num_procs = GENERATE(1, 2);
    DYNAMIC_SECTION("num_procs=" << sdpb_num_procs)
    {
      Test_Util::Test_Case_Runner runner("sdpb/run-"
                                         + std::to_string(sdpb_num_procs));
      auto args = default_args;
      run_sdpb_set_out_ck_dirs(runner, args, sdpb_num_procs);
      Test_Util::REQUIRE_Equal::diff_sdpb_output_dir(
        args["--outDir"], data_dir / "test_out_orig", 1024, 1024 / 2);
    }
  }

  SECTION("io_tests")
  {
    INFO("Check that sdpb fails on different IO errors.");
    SECTION("write_profile")
    {
      Test_Util::Test_Case_Runner runner("sdpb/io_tests/write_profile");
      auto args = default_args;
      args["--maxIterations"] = "1";
      args["--verbosity"] = "2";

      create_readonly(runner.output_dir / "ck.profiling/profiling.0");
      run_sdpb_set_out_ck_dirs(runner, args, num_procs, 1,
                               "Error when writing to");
      REQUIRE(fs::file_size(runner.output_dir / "ck.profiling/profiling.1")
              > 0);
    }

    for(std::string name :
        {"out.txt", "x_0.txt", "y.txt", "X_matrix_0.txt", "Y_matrix_0.txt"})
      {
        DYNAMIC_SECTION(name)
        {
          Test_Util::Test_Case_Runner runner("sdpb/io_tests/" + name);
          auto args = default_args;
          args["--maxIterations"] = "1";
          args["--writeSolution"] = "x,y,X,Y";

          create_readonly(runner.output_dir / "out" / name);
          std::string expected_error
            = name == "out.txt" ? "Cannot write to" : "Error when writing to";

          run_sdpb_set_out_ck_dirs(runner, args, num_procs, 1, expected_error);
        }
      }

    SECTION("input_corruption")
    {
      Test_Util::Test_Case_Runner runner("sdpb/io_tests/input_corruption");
      auto sdp_corrupted = runner.output_dir / "sdp_corrupted.zip";
      fs::create_directories(runner.output_dir);
      fs::copy(sdp_path, sdp_corrupted, fs::copy_options::recursive);
      std::ofstream os(sdp_corrupted);
      os << "any bytes to corrupt zip archive";

      auto args = default_args;
      args["--sdpDir"] = sdp_corrupted.string();
      args["--maxIterations"] = "1";

      run_sdpb_set_out_ck_dirs(runner, args, num_procs, 1);
    }

    SECTION("checkpoint_read")
    {
      Test_Util::Test_Case_Runner runner("sdpb/io_tests/checkpoint_read");
      auto args = default_args;
      args["--maxIterations"] = "1";
      args["--writeSolution"] = "x,y,X,Y";

      run_sdpb_set_out_ck_dirs(runner, args);

      INFO("Check reading from checkpoint...");
      run_sdpb_set_out_ck_dirs(runner, args);

      INFO("now use outDir as checkpoint, "
           "remove read permissions => fail to read");
      args["--checkpointDir"] = args["--outDir"];
      fs::permissions(runner.output_dir / "out" / "X_matrix_0.txt",
                      fs::perms::none);
      const Test_Util::Test_Case_Runner runner_noread
        = runner.create_nested("noread");
      runner_noread.mpi_run({"build/sdpb"}, args, num_procs, 1,
                            "Unable to open checkpoint file");
    }

    SECTION("checkpoint_corrupt")
    {
      Test_Util::Test_Case_Runner runner("sdpb/io_tests/checkpoint_corrupt");
      auto args = default_args;
      args["--maxIterations"] = "1";
      args["--writeSolution"] = "x,y,X,Y";
      run_sdpb_set_out_ck_dirs(runner, args);

      INFO("corrupt X_matrix file: remove all data after first two lines");
      auto X_matrix = (fs::path(args["--outDir"]) / "X_matrix_0.txt");
      {
        std::ifstream is(X_matrix);
        std::string line_0, line_1;
        std::getline(is, line_0);
        std::getline(is, line_1);
        is.close();

        std::ofstream os(X_matrix);
        os << line_0 << std::endl;
        os << line_1 << std::endl;
      }

      INFO("try to read checkpoint with corrupted X_matrix");
      args["--checkpointDir"] = args["--outDir"];
      const Test_Util::Test_Case_Runner runner_corrupt
        = runner.create_nested("read_corrupt");
      runner_corrupt.mpi_run({"build/sdpb"}, args, num_procs, 1,
                             "Corrupted data in file");
    }
  }
}