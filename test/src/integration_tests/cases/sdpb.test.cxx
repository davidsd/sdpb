#include "integration_tests/common.hxx"
#include <filesystem>
#include <fstream>

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
    runner.mpi_run({Test_Config::program("sdpb"), args}, num_procs,
                   required_exit_code, required_error_msg);
  };

  // create file with readonly premissions
  auto create_readonly = [](const fs::path &path) {
    if(path.has_parent_path())
      create_directories(path.parent_path());
    std::ofstream os(path);
    os << "";
    fs::permissions(path, fs::perms::others_read);
  };

#ifndef SDPB_GPU
  SECTION("gpu")
  {
    INFO("A build without the GPU backend refuses --gpu at start-up.");
    Test_Util::Test_Case_Runner runner("sdpb/gpu");
    auto args = default_args;
    args["--gpu"] = "";
    run_sdpb_set_out_ck_dirs(runner, args, 1, 1,
                             "built without the GPU backend");
  }
#endif

  SECTION("sparseB")
  {
    INFO("--sparseB stores B by its nonzeros; every result must be the same "
         "as with dense B, bit for bit.");
    // one rank, so that the comparison at full precision is deterministic
    const int np = 1;
    auto run = [&](const std::string &name,
                   const std::vector<std::string> &flags) {
      Test_Util::Test_Case_Runner runner("sdpb/sparseB/" + name);
      auto args = default_args;
      for(const auto &flag : flags)
        args[flag] = "";
      run_sdpb_set_out_ck_dirs(runner, args, np);
      return runner.output_dir / "out";
    };
    for(const bool mm : {false, true})
      DYNAMIC_SECTION((mm ? "--multimodularTrsm" : "El::Trsm"))
      {
        const std::vector<std::string> base
          = mm ? std::vector<std::string>{"--multimodularTrsm"}
               : std::vector<std::string>{};
        auto with_sparse = base;
        with_sparse.emplace_back("--sparseB");
        const std::string name = mm ? "mm" : "el";
        const auto out_dense = run(name + "/dense", base);
        const auto out_sparse = run(name + "/sparse", with_sparse);
        Test_Util::REQUIRE_Equal::diff_sdpb_output_dir(out_sparse, out_dense,
                                                       1024, 1024);
      }
    SECTION("SDP with zeros in B, 4 ranks")
    {
      INFO("dfibo test SDP (21% of B is zero), 768 bits, 4 ranks, a fixed "
           "block layout (block_timings: block 0 costs 1000, the others 1, "
           "so block 0's group has several ranks and B's blocks are "
           "distributed) and --maxSharedMemory, 40 iterations: --sparseB "
           "gives dense B's output files bit for bit.");
      const auto dfibo_sdp = Test_Config::test_data_dir / "end-to-end_tests"
                             / "dfibo-0-0-j=3-c=3.0000-d=3-s=6" / "output"
                             / "sdp";
      size_t num_blocks = 0;
      for(const auto &entry : fs::directory_iterator(dfibo_sdp))
        if(entry.path().filename().string().rfind("block_info_", 0) == 0)
          ++num_blocks;
      REQUIRE(num_blocks > 1);
      auto run_dfibo = [&](const std::string &name,
                           const std::vector<std::string> &flags) {
        Test_Util::Test_Case_Runner runner("sdpb/sparseB/dfibo/" + name);
        // the layout from a file, not from a timing run
        create_directories(runner.output_dir / "ck");
        {
          std::ofstream os(runner.output_dir / "ck" / "block_timings");
          for(size_t b = 0; b < num_blocks; ++b)
            os << (b == 0 ? 1000 : 1) << "\n";
        }
        Test_Util::Test_Case_Runner::Named_Args_Map args{
          {"--sdpDir", dfibo_sdp.string()},
          {"--precision", "768"},
          {"--noFinalCheckpoint", ""},
          {"--maxIterations", "40"},
          {"--maxSharedMemory", "1G"},
          {"--initialMatrixScalePrimal", "1e10"},
          {"--initialMatrixScaleDual", "1e10"},
          {"--writeSolution", "x,y"}};
        for(const auto &flag : flags)
          args[flag] = "";
        run_sdpb_set_out_ck_dirs(runner, args, 4);
        return runner.output_dir / "out";
      };
      for(const bool mm : {false, true})
        DYNAMIC_SECTION((mm ? "--multimodularTrsm" : "El::Trsm"))
        {
          const std::vector<std::string> base
            = mm ? std::vector<std::string>{"--multimodularTrsm"}
                 : std::vector<std::string>{};
          auto with_sparse = base;
          with_sparse.emplace_back("--sparseB");
          const std::string name = mm ? "mm" : "el";
          const auto out_dense = run_dfibo(name + "/dense", base);
          const auto out_sparse = run_dfibo(name + "/sparse", with_sparse);
          Test_Util::REQUIRE_Equal::diff_sdpb_output_dir(out_sparse,
                                                         out_dense, 768, 768);
        }
    }
    SECTION("SDPA input")
    {
      Test_Util::Test_Case_Runner runner("sdpb/sparseB/sdpa");
      auto args = default_args;
      args["--sdpDir"] = (Test_Config::test_data_dir / "end-to-end_tests"
                          / "sdpa" / "example1" / "input" / "example1.dat-s")
                           .string();
      args["--precision"] = "200";
      args["--sparseB"] = "";
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
        runner_extra_ck_file.mpi_run({Test_Config::program("sdpb"), args},
                                     num_procs, 1,
                                     "unexpected checkpoint file");
      }

      {
        INFO("Remove checkpoint file for rank=1 => loading error");
        fs::remove(args["--checkpointDir"] + "/checkpoint_2_1");
        const Test_Util::Test_Case_Runner runner_missing_ck_file
          = runner.create_nested("missing_ck_file");
        runner_missing_ck_file.mpi_run({Test_Config::program("sdpb"), args},
                                       num_procs, 1,
                                       "Cannot find checkpoint file");
      }

      INFO("now use outDir as checkpoint, "
           "remove read permissions => fail to read");
      args["--checkpointDir"] = args["--outDir"];
      fs::permissions(runner.output_dir / "out" / "X_matrix_0.txt",
                      fs::perms::none);
      const Test_Util::Test_Case_Runner runner_noread
        = runner.create_nested("noread");
      runner_noread.mpi_run({Test_Config::program("sdpb"), args}, num_procs, 1,
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
      runner_corrupt.mpi_run({Test_Config::program("sdpb"), args}, num_procs,
                             1, "Corrupted data in file");
    }
  }
}