#include "integration_tests/common.hxx"

// Realistic end-to-end test for sdp2input + sdpb
// JSON input taken from  "SingletScalar_cT_test_nmax6" and
// "SingletScalarAllowed_test_nmax6"
// https://gitlab.com/davidsd/scalars-3d/-/blob/master/src/Projects/Scalars3d/SingletScalar2020.hs
// Test data is generated with SDPB 2.5.1 on Caltech cluster.
// Note that on different machines results can vary due to rounding errors,
// depending on GMP/MPFR version etc.

namespace
{
  void end_to_end_test(const std::string &name,
                       const std::string &default_sdpb_args, int num_procs,
                       int precision, int sdp2input_final_precision,
                       int sdpb_final_precision,
                       const std::vector<std::string> &out_txt_keys = {})
  {
    Test_Util::Test_Case_Runner runner(name);
    const auto &data_dir = runner.data_dir;
    const auto &output_dir = runner.output_dir;

    auto sdp_orig_zip = (data_dir / ("sdp.orig.zip")).string();
    auto sdp_zip = (output_dir / "sdp.zip").string();

    // sdp2input
    {
      Test_Util::Test_Case_Runner::Named_Args_Map args{
        {"--input", (data_dir / "json" / "file_list.nsv").string()},
        {"--output", sdp_zip},
        {"--precision", std::to_string(precision)}};

      runner.create_nested("sdp2input/run")
        .mpi_run({"build/sdp2input"}, args, num_procs);

      // sdp2input runs with --precision=<precision>
      // We check test output up to lower precision=<sdp2input_final_precision>
      // in order to neglect unimportant rounding errors
      Test_Util::REQUIRE_Equal::diff_sdp_zip(
        sdp_zip, sdp_orig_zip, sdp2input_final_precision,
        runner.create_nested("sdp2input/diff"));
    }

    // sdpb
    {
      Test_Util::Test_Case_Runner::Named_Args_Map args{
        {"--procsPerNode", std::to_string(num_procs)},
        {"--precision", std::to_string(precision)},
        {"--sdpDir", sdp_zip},
        {"--outDir", (output_dir / "out").string()},
        {"--checkpointDir", (output_dir / "ck").string()}};
      runner.create_nested("sdpb/run")
        .mpi_run({"build/sdpb", default_sdpb_args}, args, num_procs);

      // SDPB runs with --precision=<precision>
      // We check test output up to lower precision=<sdpb_final_precision>
      // in order to neglect unimportant rounding errors
      Test_Util::REQUIRE_Equal::diff_sdpb_output_dir(
        output_dir / "out", data_dir / "out", sdpb_final_precision, {},
        out_txt_keys);
    }
  }
}

TEST_CASE("end-to-end_tests")
{
  INFO("End-to-end tests for sdp2input + sdpb");
  INFO("Test data is generated with SDPB 2.5.1 on Caltech cluster.");
  INFO("On different machines results can vary due to rounding errors, "
       "depending on GMP/MPFR version etc");
  int num_procs = 6;
  int precision = 768;
  int sdp2input_final_precision = 720;
  std::string name = "end-to-end_tests";

  SECTION("SingletScalar_cT_test_nmax6/primal_dual_optimal")
  {
    INFO("SingletScalar_cT_test_nmax6 from "
         "https://gitlab.com/davidsd/scalars-3d/-/blob/master/src/Projects/"
         "Scalars3d/SingletScalar2020.hs");
    INFO("SDPB should find primal-dual optimal solution.");
    name += "/SingletScalar_cT_test_nmax6/primal_dual_optimal";
    std::string default_sdpb_args
      = "--checkpointInterval 3600 --maxRuntime 1340 "
        "--dualityGapThreshold 1.0e-30 --primalErrorThreshold 1.0e-30 "
        "--dualErrorThreshold 1.0e-30 --initialMatrixScalePrimal 1.0e20 "
        "--initialMatrixScaleDual 1.0e20 --feasibleCenteringParameter 0.1 "
        "--infeasibleCenteringParameter 0.3 --stepLengthReduction 0.7 "
        "--maxComplementarity 1.0e100 --maxIterations 1000 --verbosity 1 "
        "--procGranularity 1 --writeSolution y";
    int sdpb_final_precision = 600;
    end_to_end_test(name, default_sdpb_args, num_procs, precision,
                    sdp2input_final_precision, sdpb_final_precision);
  }

  SECTION("SingletScalarAllowed_test_nmax6")
  {
    INFO("SingletScalarAllowed_test_nmax6 from "
         "https://gitlab.com/davidsd/scalars-3d/-/blob/master/src/Projects/"
         "Scalars3d/SingletScalar2020.hs");
    name += "/SingletScalarAllowed_test_nmax6";
    std::string default_sdpb_args
      = "--checkpointInterval 3600 --maxRuntime 1341 "
        "--dualityGapThreshold 1.0e-30 --primalErrorThreshold 1.0e-200 "
        "--dualErrorThreshold 1.0e-200 --initialMatrixScalePrimal 1.0e20 "
        "--initialMatrixScaleDual 1.0e20 --feasibleCenteringParameter 0.1 "
        "--infeasibleCenteringParameter 0.3 --stepLengthReduction 0.7 "
        "--maxComplementarity 1.0e100 --maxIterations 1000 --verbosity 1 "
        "--procGranularity 1 --writeSolution y "
        "--detectPrimalFeasibleJump --detectDualFeasibleJump";
    int sdpb_final_precision = 610;

    SECTION("primal_feasible_jump")
    {
      INFO("SDPB should detect primal feasible jump.");
      name += "/primal_feasible_jump";
      std::vector<std::string> out_txt_keys
        = {"terminateReason", "primalObjective", "dualObjective", "dualityGap",
           "dualError"};
      end_to_end_test(name, default_sdpb_args, num_procs, precision,
                      sdp2input_final_precision, sdpb_final_precision,
                      out_txt_keys);
    }
    SECTION("dual_feasible_jump")
    {
      INFO("SDPB should detect dual feasible jump.");
      name += "/dual_feasible_jump";
      std::vector<std::string> out_txt_keys
        = {"terminateReason", "primalObjective", "dualObjective", "dualityGap",
           "primalError"};
      end_to_end_test(name, default_sdpb_args, num_procs, precision,
                      sdp2input_final_precision, sdpb_final_precision,
                      out_txt_keys);
    }
  }
}
