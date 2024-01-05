#include "integration_tests/common.hxx"

#include <filesystem>

using namespace std::string_literals;
namespace fs = std::filesystem;

TEST_CASE("outer_limits")
{
  INFO("Simple outer_limits test");
  Test_Util::Test_Case_Runner outer_runner("outer_limits");
  fs::path data_dir = outer_runner.data_dir; //test/data/outer_limits

  int num_procs = GENERATE(1, 2, 6);
  DYNAMIC_SECTION("num_procs=" << num_procs)
  {
    Test_Util::Test_Case_Runner runner
      = outer_runner.create_nested("mpirun-" + std::to_string(num_procs));

    fs::path output_dir = runner.output_dir;

    Test_Util::Test_Case_Runner::Named_Args_Map args{
      {"--functions", (data_dir / "toy_functions.json").string()},
      {"--out", (output_dir / "toy_functions_out.json").string()},
      {"--checkpointDir", (output_dir / "ck").string()},
      {"--points", (data_dir / "toy_functions_points.json").string()},
      {"--precision", "128"},
      {"--dualityGapThreshold", "1e-10"},
      {"--primalErrorThreshold", "1e-10"},
      {"--dualErrorThreshold", "1e-10"},
      {"--initialMatrixScalePrimal", "1e1"},
      {"--initialMatrixScaleDual", "1e1"},
      {"--maxIterations", "1000"},
      {"--verbosity", "1"},
    };

    runner.mpi_run({"build/outer_limits"}, args, num_procs);

    auto out = output_dir / "toy_functions_out.json";
    auto out_orig = data_dir / "toy_functions_out_orig.json";
    Test_Util::REQUIRE_Equal::diff_outer_limits(out, out_orig, 128, 128);
  }
}