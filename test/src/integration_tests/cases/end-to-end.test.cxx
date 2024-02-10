#include "integration_tests/common.hxx"

// Realistic end-to-end test for pmp2sdp + sdpb
// JSON input taken from  "SingletScalar_cT_test_nmax6" and
// "SingletScalarAllowed_test_nmax6"
// https://gitlab.com/davidsd/scalars-3d/-/blob/master/src/Projects/Scalars3d/SingletScalar2020.hs
// Test data is generated with SDPB 2.5.1 on Caltech cluster.
// Note that on different machines results can vary due to rounding errors,
// depending on GMP/MPFR version etc.

using namespace Test_Util;
using namespace Test_Util::REQUIRE_Equal;
using Named_Args_Map = Test_Case_Runner::Named_Args_Map;

namespace fs = std::filesystem;

namespace
{
  void end_to_end_test(const std::string &name, int num_procs, int precision,
                       const std::string &default_sdpb_args = "",
                       const Named_Args_Map &pmp2sdp_args = {},
                       const std::vector<std::string> &out_txt_keys = {},
                       bool check_sdp_normalization = true)
  {
    int diff_precision = precision / 2;

    const auto data_dir
      = Test_Config::test_data_dir / "end-to-end_tests" / name;
    auto data_input_dir = data_dir / "input";
    auto data_output_dir = data_dir / "output";

    std::string sdp_format;
    if(pmp2sdp_args.find("--outputFormat") != pmp2sdp_args.end())
      sdp_format = pmp2sdp_args.at("--outputFormat");
    auto format_description = sdp_format.empty() ? "default(bin)" : sdp_format;
    CAPTURE(sdp_format);
    CAPTURE(format_description);

    std::vector<fs::path> pmp_paths;
    for(const auto &it : fs::directory_iterator(data_input_dir))
      {
        if(is_regular_file(it.path()))
          pmp_paths.push_back(it.path());
      }
    REQUIRE(!pmp_paths.empty());

    for(const auto &pmp_path : pmp_paths)
      DYNAMIC_SECTION(pmp_path.filename().string())
      {
        std::string runner_name = name + "/" + pmp_path.filename().string();
        if(!sdp_format.empty())
          runner_name += "/format=" + sdp_format;
        Test_Case_Runner runner(runner_name);
        const auto &output_dir = runner.output_dir;

        bool zip = pmp2sdp_args.find("--zip") != pmp2sdp_args.end();
        auto sdp_path = (output_dir / (zip ? "sdp.zip" : "sdp")).string();
        // pmp2sdp
        {
          INFO("run pmp2sdp");
          auto args = pmp2sdp_args;

          args.insert({{"--input", pmp_path.string()},
                       {"--output", sdp_path},
                       {"--precision", std::to_string(precision)}});

          runner.create_nested("pmp2sdp").mpi_run({"build/pmp2sdp"}, args,
                                                  num_procs);

          // pmp2sdp runs with --precision=<precision>
          // We check test output up to lower precision=<diff_precision>
          // in order to neglect unimportant rounding errors
          auto sdp_orig = data_output_dir / "sdp";
          diff_sdp(sdp_path, sdp_orig, precision, diff_precision,
                   runner.create_nested("sdp.diff"), check_sdp_normalization);
        }

        // sdpb
        {
          Named_Args_Map args{
            {"--precision", std::to_string(precision)},
            {"--sdpDir", sdp_path},
            {"--outDir", (output_dir / "out").string()},
            {"--checkpointDir", (output_dir / "ck").string()}};
          runner.create_nested("sdpb").mpi_run(
            {"build/sdpb", default_sdpb_args}, args, num_procs);

          // SDPB runs with --precision=<precision>
          // We check test output up to lower precision=<sdpb_output_diff_precision>
          // in order to neglect unimportant rounding errors
          diff_sdpb_output_dir(output_dir / "out", data_output_dir / "out",
                               precision, diff_precision, {}, out_txt_keys);
        }

        if(exists(data_output_dir / "spectrum.json"))
          {
            Named_Args_Map args{
              {"--input", pmp_path.string()},
              {"--solution", (output_dir / "out").string()},
              {"--threshold", "1e-10"},
              {"--output", (output_dir / "spectrum.json").string()},
              {"--precision", std::to_string(precision)},
            };
            runner.create_nested("spectrum")
              .mpi_run({"build/spectrum"}, args, num_procs);
            diff_spectrum(output_dir / "spectrum.json",
                          data_output_dir / "spectrum.json", precision,
                          diff_precision);
          }
      }
  }

  Named_Args_Map
  build_pmp2sdp_args(const std::string &sdp_format, const bool zip = false)
  {
    Named_Args_Map result;
    if(!sdp_format.empty())
      result["--outputFormat"] = sdp_format;
    if(zip)
      result["--zip"] = "";
    return result;
  }
}

TEST_CASE("end-to-end_tests")
{
  INFO("End-to-end tests for pmp2sdp + sdpb");
  INFO("On different machines results can vary due to rounding errors, "
       "depending on GMP/MPFR version etc");
  int num_procs = 6;
  int precision = 768;

  SECTION("1d")
  {
    INFO("SDPB test for a simple one-dimensional problem from SDPB Manual:");
    INFO("maximize (-y) s.t. (1 + x^4 + y * (x^4 / 12 + x^2)) >= 0)");

    num_procs = 2;
    // Do not check normalization.json because it is absent for XML
    bool check_sdp_normalization = false;
    end_to_end_test("1d", num_procs, precision, {}, {}, {},
                    check_sdp_normalization);
  }

  SECTION("dfibo-0-0-j=3-c=3.0000-d=3-s=6")
  {
    INFO("pmp2sdp+sdpb test for https://github.com/davidsd/sdpb/issues/124");
    INFO("sdp contains block with empty bilinear_bases_odd, "
         "which caused a bug.");
    INFO("Test data from Harvard cluster, gmp/6.2.1 mpfr/4.2.0");
    std::string default_sdpb_args
      = "--findDualFeasible --findPrimalFeasible "
        "--initialMatrixScalePrimal 1e10 --initialMatrixScaleDual 1e10 "
        "--maxComplementarity 1e30 --dualErrorThreshold 1e-10 "
        "--primalErrorThreshold 1e-153 --maxRuntime 259200 "
        "--checkpointInterval 3600 --maxIterations 1000 "
        "--feasibleCenteringParameter=0.1 --infeasibleCenteringParameter=0.3 "
        "--stepLengthReduction=0.7";
    for(std::string sdp_format : {"", "bin", "json"})
      {
        DYNAMIC_SECTION(
          "format=" << (sdp_format.empty() ? "default(bin)" : sdp_format))
        {
          // write sdp to zip instead of plain directory
          bool zip = true;
          end_to_end_test("dfibo-0-0-j=3-c=3.0000-d=3-s=6", num_procs,
                          precision, default_sdpb_args,
                          build_pmp2sdp_args(sdp_format, zip));
        }
      }
  }

  SECTION("SingletScalar_cT_test_nmax6/primal_dual_optimal")
  {
    INFO("SingletScalar_cT_test_nmax6 from "
         "https://gitlab.com/davidsd/scalars-3d/-/blob/master/src/Projects/"
         "Scalars3d/SingletScalar2020.hs");
    INFO("Test data is generated with SDPB 2.5.1 on Caltech cluster.");
    INFO("SDPB should find primal-dual optimal solution.");
    const auto name = "SingletScalar_cT_test_nmax6/primal_dual_optimal";
    std::string default_sdpb_args
      = "--checkpointInterval 3600 --maxRuntime 1340 "
        "--dualityGapThreshold 1.0e-30 --primalErrorThreshold 1.0e-30 "
        "--dualErrorThreshold 1.0e-30 --initialMatrixScalePrimal 1.0e20 "
        "--initialMatrixScaleDual 1.0e20 --feasibleCenteringParameter 0.1 "
        "--infeasibleCenteringParameter 0.3 --stepLengthReduction 0.7 "
        "--maxComplementarity 1.0e100 --maxIterations 1000 --verbosity 1 "
        "--procGranularity 1 --writeSolution x,y,z";
    // This test is slow, we don't want to run it twice
    // json/bin correctness is checked by other tests,
    // so we use only binary SDP here
    end_to_end_test(name, num_procs, precision, default_sdpb_args);
  }

  SECTION("SingletScalarAllowed_test_nmax6")
  {
    INFO("SingletScalarAllowed_test_nmax6 from "
         "https://gitlab.com/davidsd/scalars-3d/-/blob/master/src/Projects/"
         "Scalars3d/SingletScalar2020.hs");
    INFO("Test data is generated with SDPB 2.5.1 on Caltech cluster.");
    std::string name = "SingletScalarAllowed_test_nmax6";
    std::string default_sdpb_args
      = "--checkpointInterval 3600 --maxRuntime 1341 "
        "--dualityGapThreshold 1.0e-30 --primalErrorThreshold 1.0e-200 "
        "--dualErrorThreshold 1.0e-200 --initialMatrixScalePrimal 1.0e20 "
        "--initialMatrixScaleDual 1.0e20 --feasibleCenteringParameter 0.1 "
        "--infeasibleCenteringParameter 0.3 --stepLengthReduction 0.7 "
        "--maxComplementarity 1.0e100 --maxIterations 1000 --verbosity 1 "
        "--procGranularity 1 --writeSolution y,z "
        "--detectPrimalFeasibleJump --detectDualFeasibleJump";

    SECTION("primal_feasible_jump")
    {
      INFO("SDPB should detect primal feasible jump.");
      name += "/primal_feasible_jump";
      std::vector<std::string> out_txt_keys
        = {"terminateReason", "primalObjective", "dualObjective", "dualityGap",
           "dualError"};
      end_to_end_test(name, num_procs, precision, default_sdpb_args, {},
                      out_txt_keys);
    }
    SECTION("dual_feasible_jump")
    {
      INFO("SDPB should detect dual feasible jump.");
      name += "/dual_feasible_jump";
      std::vector<std::string> out_txt_keys
        = {"terminateReason", "primalObjective", "dualObjective", "dualityGap",
           "primalError"};

      end_to_end_test(name, num_procs, precision, default_sdpb_args, {},
                      out_txt_keys);
    }
  }
}
