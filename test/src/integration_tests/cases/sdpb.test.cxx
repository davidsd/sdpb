#include "integration_tests/common.hxx"
#include "sdpb_util/Boost_Float.hxx"
#include <rapidjson/document.h>
#include <rapidjson/istreamwrapper.h>
#include <filesystem>
#include <fstream>
#include <iterator>

namespace fs = std::filesystem;
using Test_Util::REQUIRE_Equal::diff;

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
                // At --verbosity debug the profile keeps the timers that
                // have finished, not only the ones still running
                std::ifstream is(profiling_path);
                const std::string profile(
                  (std::istreambuf_iterator<char>(is)), {});
                CAPTURE(profiling_path);
                REQUIRE(profile.find("\"sdpb.solve.read_sdp\"")
                        != std::string::npos);
                REQUIRE(profile.find("\"sdpb.solve.run.iter_1\"")
                        != std::string::npos);
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
namespace
{
  // Entries of a vector written by write_matrix() ("<N> 1", then N entries)
  std::vector<std::string> read_vector_txt(const fs::path &path)
  {
    CAPTURE(path);
    std::ifstream is(path);
    size_t height = 0, width = 0;
    is >> height >> width;
    REQUIRE(is.good());
    REQUIRE(width == 1);
    std::vector<std::string> entries(height);
    for(auto &entry : entries)
      is >> entry;
    REQUIRE(!is.fail());
    return entries;
  }

  // "normalization" of normalization.json, as the strings in the file
  std::vector<std::string> read_normalization_json(const fs::path &path)
  {
    CAPTURE(path);
    std::ifstream is(path);
    REQUIRE(is.good());
    rapidjson::IStreamWrapper wrapper(is);
    rapidjson::Document document;
    document.ParseStream(wrapper);
    std::vector<std::string> result;
    for(const auto &value : document["normalization"].GetArray())
      result.emplace_back(value.GetString());
    return result;
  }
}

TEST_CASE("sdpb_z_normalization")
{
  INFO("z.txt must satisfy n.z == 1 as written: z_k = (1 - sum_{i!=k} "
       "n_i z_i) / n_k formed and written at SDPB's precision gives n.z far "
       "from 1 when |n_i z_i| >> 2^precision, e.g. after a dual feasible "
       "jump.");
  INFO("Here: the 1d problem with a normalization of many digits, "
       "--precision=64, and a hand-made y = -1.23e60 loaded from a text "
       "checkpoint; SDPB stops at once (--maxIterations=0) and writes z.");
  Test_Util::Test_Case_Runner runner("sdpb_z_normalization");
  const auto &dir = runner.output_dir;
  fs::create_directories(dir);
  const int num_procs = 2;

  // The 1d problem with normalization 1/sqrt(2), 1/pi (200 digits; k = 0)
  // and a huge y (input/y.txt)
  const auto input_dir = Test_Config::test_data_dir / "end-to-end_tests"
                         / "1d-z-normalization" / "input";
  const auto pmp_path = input_dir / "pmp.json";
  const auto sdp_path = dir / "sdp";
  runner.create_nested("pmp2sdp").mpi_run(
    {Test_Config::program("pmp2sdp"),
     {{"--input", pmp_path.string()},
      {"--output", sdp_path.string()},
      {"--precision", "1024"},
      // json, so that SDPB can read it at another precision
      {"--outputFormat", "json"}}},
    num_procs);

  const Test_Util::Test_Case_Runner::Named_Args_Map args{
    {"--sdpDir", sdp_path.string()},
    {"--precision", "64"},
    {"--noFinalCheckpoint", ""},
  };
  // One iteration, to get a text checkpoint (x, y, X, Y)
  const auto out_1 = dir / "out_1";
  {
    auto a = args;
    a["--maxIterations"] = "1";
    a["--writeSolution"] = "x,y,X,Y";
    a["--outDir"] = out_1.string();
    a["--checkpointDir"] = (dir / "ck_1").string();
    runner.create_nested("sdpb_1").mpi_run({Test_Config::program("sdpb"), a},
                                          num_procs);
  }
  // Replace y by a huge one
  fs::copy_file(input_dir / "y.txt", out_1 / "y.txt",
                fs::copy_options::overwrite_existing);
  const auto out_2 = dir / "out_2";
  {
    auto a = args;
    a["--maxIterations"] = "0";
    a["--writeSolution"] = "y,z";
    a["--outDir"] = out_2.string();
    a["--checkpointDir"] = out_1.string();
    const auto sdpb_2 = runner.create_nested("sdpb_2");
    sdpb_2.mpi_run({Test_Config::program("sdpb"), a}, num_procs);

    INFO("the z.txt self-check prints no warning");
    std::ifstream is(sdpb_2.stderr_path);
    const std::string stderr_text((std::istreambuf_iterator<char>(is)),
                                  std::istreambuf_iterator<char>());
    CAPTURE(stderr_text);
    REQUIRE(stderr_text.find("n.z") == std::string::npos);
  }

  const auto z = read_vector_txt(out_2 / "z.txt");
  const auto y = read_vector_txt(out_2 / "y.txt");
  const auto n = read_normalization_json(sdp_path / "normalization.json");
  REQUIRE(z.size() == 2);
  REQUIRE(y.size() == 1);
  REQUIRE(n.size() == 2);
  CAPTURE(z, y, n);

  INFO("z_i (i != k = 0) are written as y.txt writes y");
  DIFF(std::vector<std::string>(z.begin() + 1, z.end()), y);

  INFO("|z_1| must be the huge y");
  Scoped_Boost_Float_Precision precision(5000);
  REQUIRE(abs(Boost_Float(z[1])) > 1e59);

  INFO("n.z == 1 for the written decimals, to 2^-64 (at 5000 bits)");
  Boost_Float nz = 0;
  for(size_t i = 0; i < z.size(); ++i)
    nz += Boost_Float(n[i]) * Boost_Float(z[i]);
  const Boost_Float error = abs(nz - 1);
  CAPTURE(error);
  REQUIRE(error <= ldexp(Boost_Float(1), -64));
}
