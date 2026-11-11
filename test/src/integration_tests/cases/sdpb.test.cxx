#include "integration_tests/common.hxx"
#include "integration_tests/util/Float.hxx"
#include <filesystem>
#include <fstream>
#include <rapidjson/istreamwrapper.h>
#include <rapidjson/document.h>

namespace fs = std::filesystem;

TEST_CASE("sdpb")
{
  INFO("SDPB tests for a simple one-dimensional problem:");
  INFO("maximize (-y) s.t. (1 + x^4 + y * (x^4 / 12 + x^2)) >= 0)");
  // default sdp input file
  auto sdp_path = Test_Config::test_data_dir / "end-to-end_tests" / "1d"
                  / "output" / "sdp";
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
    runner.mpi_run({"build/sdpb", args}, num_procs, required_exit_code,
                   required_error_msg);
  };

  // create file with readonly premissions
  auto create_readonly = [](const fs::path &path) {
    if(path.has_parent_path())
      create_directories(path.parent_path());
    std::ofstream os(path);
    os << "";
    fs::permissions(path, fs::perms::others_read);
  };

  SECTION("reconstructXYmu")
  {
    INFO("A text checkpoint with only x and y: refused without "
         "--reconstructXYmu; with it, X = sum_p A_p x_p, Y = mu X^{-1}, and "
         "the run reaches the same answer as a run from the start.");
    const int np = 1;
    // x and y of the first primal feasible iterate, as text. There
    // X = sum_p A_p x_p (to the primal error) is well inside the cone; an
    // earlier, primal infeasible x can give an indefinite sum.
    Test_Util::Test_Case_Runner runner_xy("sdpb/reconstructXYmu/xy");
    auto args_xy = default_args;
    args_xy["--findPrimalFeasible"] = "";
    args_xy["--writeSolution"] = "x,y";
    run_sdpb_set_out_ck_dirs(runner_xy, args_xy, np);
    const auto xy_dir = runner_xy.output_dir / "out";
    REQUIRE(fs::exists(xy_dir / "x_0.txt"));
    REQUIRE(!fs::exists(xy_dir / "X_matrix_0.txt"));

    Test_Util::Test_Case_Runner runner_full("sdpb/reconstructXYmu/full");
    auto args_full = default_args;
    run_sdpb_set_out_ck_dirs(runner_full, args_full, np);

    {
      INFO("without the option the checkpoint is refused");
      Test_Util::Test_Case_Runner runner("sdpb/reconstructXYmu/refused");
      auto args = default_args;
      args["--initialCheckpointDir"] = xy_dir.string();
      run_sdpb_set_out_ck_dirs(runner, args, np, 1,
                               "has no X_matrix/Y_matrix files");
    }
    {
      Test_Util::Test_Case_Runner runner("sdpb/reconstructXYmu/restart");
      auto args = default_args;
      args["--initialCheckpointDir"] = xy_dir.string();
      // about the mu of that iterate
      args["--reconstructXYmu"] = "1e15";
      run_sdpb_set_out_ck_dirs(runner, args, np);
      Test_Util::REQUIRE_Equal::diff_sdpb_output_dir(
        runner.output_dir / "out", runner_full.output_dir / "out", 1024, 64,
        {"out.txt"}, {"terminateReason", "primalObjective", "dualObjective"});
      {
        INFO("The first iterate is the reconstructed one: X = sum_p A_p x_p "
             "(primal residue P-err ~ 0) and mu = Tr(XY)/dim = 1e15.");
        std::ifstream is(runner.output_dir / "out" / "iterations.json");
        rapidjson::IStreamWrapper wrapper(is);
        rapidjson::Document doc;
        doc.ParseStream(wrapper);
        REQUIRE(doc.IsArray());
        REQUIRE(doc.Size() > 0);
        const auto &first = doc[0];
        Float_Binary_Precision prec(1024, 1024);
        const El::BigFloat mu(first["mu"].GetString());
        const El::BigFloat P_err(first["P-err"].GetString());
        CAPTURE(mu);
        CAPTURE(P_err);
        REQUIRE(El::Abs(mu / El::BigFloat("1e15") - 1)
                < El::BigFloat("1e-200"));
        REQUIRE(P_err < El::BigFloat("1e-200"));
      }
    }
    {
      INFO("x, y and X/Y for some blocks only: refused");
      Test_Util::Test_Case_Runner runner_part("sdpb/reconstructXYmu/partial");
      auto args_part = default_args;
      args_part["--maxIterations"] = "2";
      args_part["--writeSolution"] = "x,y,X,Y";
      run_sdpb_set_out_ck_dirs(runner_part, args_part, np);
      const auto part_dir = runner_part.output_dir / "out";
      REQUIRE(fs::exists(part_dir / "Y_matrix_0.txt"));
      fs::remove(part_dir / "Y_matrix_0.txt");
      Test_Util::Test_Case_Runner runner("sdpb/reconstructXYmu/partial/run");
      auto args = default_args;
      args["--initialCheckpointDir"] = part_dir.string();
      run_sdpb_set_out_ck_dirs(runner, args, np, 1,
                               "for some blocks but not for others");
    }
    {
      INFO("SDPA input is refused");
      Test_Util::Test_Case_Runner runner("sdpb/reconstructXYmu/sdpa");
      auto args = default_args;
      args["--sdpDir"] = (Test_Config::test_data_dir / "end-to-end_tests"
                          / "sdpa" / "example1" / "input" / "example1.dat-s")
                           .string();
      args["--precision"] = "200";
      args["--reconstructXYmu"] = "1e-10";
      run_sdpb_set_out_ck_dirs(runner, args, np, 1, "not for SDPA input");
    }
  }

  SECTION("io_tests")
  {
    INFO("Check that sdpb fails on different IO errors.");
    SECTION("write_profile")
    {
      int max_iterations = GENERATE(1, 2);
      DYNAMIC_SECTION("--maxIterations=" << max_iterations)
      {
        Test_Util::Test_Case_Runner runner(
          "sdpb/io_tests/write_profile/--maxIterations="
          + std::to_string(max_iterations));
        auto args = default_args;
        args["--maxIterations"] = std::to_string(max_iterations);
        args["--verbosity"] = "2";

        create_readonly(runner.output_dir / "ck.profiling/profiling.0");
        run_sdpb_set_out_ck_dirs(runner, args, num_procs);

        {
          INFO(
            "ck.profiling/ exists, SDPB should rename it to ck.profiling.0/");
          INFO("Then SDPB should create "
               "ck.profiling.1/profiling.* for timing run and "
               "ck.profiling/profiling.* for actual run.");
          for(const std::string dir_suffix : {"", ".1"})
            for(const size_t rank : {0, 1})
              {
                auto profiling_path = runner.output_dir
                                      / ("ck.profiling" + dir_suffix)
                                      / ("profiling." + std::to_string(rank));
                REQUIRE(fs::file_size(profiling_path) > 0);
              }
        }
        {
          REQUIRE(fs::file_size(runner.output_dir / "ck/block_timings") > 0);
        }
      }
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
      // We'll need a checkpoint in this case
      args.erase("--noFinalCheckpoint");

      run_sdpb_set_out_ck_dirs(runner, args);

      INFO("Check reading from checkpoint...");
      run_sdpb_set_out_ck_dirs(runner, args);

      {
        INFO("Add extra checkpoint file => loading error");
        std::ofstream os(args["--checkpointDir"] + "/checkpoint_2_2");
        os << "any bytes to checkpoint_2_2";
        const Test_Util::Test_Case_Runner runner_extra_ck_file
          = runner.create_nested("extra_ck_file");
        runner_extra_ck_file.mpi_run({"build/sdpb", args}, num_procs, 1,
                                     "unexpected checkpoint file");
      }

      {
        INFO("Remove checkpoint file for rank=1 => loading error");
        fs::remove(args["--checkpointDir"] + "/checkpoint_2_1");
        const Test_Util::Test_Case_Runner runner_missing_ck_file
          = runner.create_nested("missing_ck_file");
        runner_missing_ck_file.mpi_run({"build/sdpb", args}, num_procs, 1,
                                       "Cannot find checkpoint file");
      }

      INFO("now use outDir as checkpoint, "
           "remove read permissions => fail to read");
      args["--checkpointDir"] = args["--outDir"];
      fs::permissions(runner.output_dir / "out" / "X_matrix_0.txt",
                      fs::perms::none);
      const Test_Util::Test_Case_Runner runner_noread
        = runner.create_nested("noread");
      runner_noread.mpi_run({"build/sdpb", args}, num_procs, 1,
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
      runner_corrupt.mpi_run({"build/sdpb", args}, num_procs, 1,
                             "Corrupted data in file");
    }
  }
}