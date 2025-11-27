#include "integration_tests/common.hxx"

#include <boost/program_options/parsers.hpp>

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
  struct End_To_End_Test
  {
    std::string name;
    int num_procs = 6;
    int precision = 768;
    Named_Args_Map pmp2sdp_args;
    std::vector<std::string> default_sdpb_args;
    std::vector<std::string> sdpb_out_filenames;
    std::vector<std::string> sdpb_out_txt_keys;
    bool check_sdp = true;
    bool check_sdp_normalization = true;
    bool run_sdpb_twice = false;

    explicit End_To_End_Test(std::string name) : name(std::move(name)) {}

    void run()
    {
      int diff_precision = precision / 2;

      const auto data_dir
        = Test_Config::test_data_dir / "end-to-end_tests" / name;
      auto data_input_dir = data_dir / "input";
      auto data_output_dir = data_dir / "output";

      for(const auto &it : fs::directory_iterator(data_output_dir))
        {
          INFO("Check filenames in data_output_dir");
          CAPTURE(data_output_dir);
          auto filename = it.path().filename();
          CAPTURE(filename);
          if(!(filename == "sdp" || filename == "out"
               || filename == "spectrum.json" || filename == "iterations.json"
               || filename == "iterations.0.json"
               || filename == "iterations.1.json"))
            {
              FAIL("Unexpected file: " << it.path());
            }
        }

      std::string sdp_format;
      if(pmp2sdp_args.find("--outputFormat") != pmp2sdp_args.end())
        sdp_format = pmp2sdp_args.at("--outputFormat");
      auto format_description
        = sdp_format.empty() ? "default(bin)" : sdp_format;
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
          std::string runner_name
            = "end-to-end_tests/" + name + "/" + pmp_path.filename().string();
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

            runner.create_nested("pmp2sdp").mpi_run({"build/pmp2sdp", args},
                                                    num_procs);

            if(check_sdp)
              {
                // pmp2sdp runs with --precision=<precision>
                // We check test output up to lower precision=<diff_precision>
                // in order to neglect unimportant rounding errors
                auto sdp_orig = data_output_dir / "sdp";
                diff_sdp(sdp_path, sdp_orig, precision, diff_precision,
                         runner.create_nested("sdp.diff"),
                         check_sdp_normalization);
              }
          }

          // sdpb
          {
            Named_Args_Map args{
              {"--precision", std::to_string(precision)},
              {"--sdpDir", sdp_path},
              {"--outDir", (output_dir / "out").string()},
              {"--checkpointDir", (output_dir / "ck").string()}};
            runner.create_nested("sdpb").mpi_run(
              {"build/sdpb", default_sdpb_args, args}, num_procs);
            if(run_sdpb_twice)
              {
                runner.create_nested("sdpb-2").mpi_run(
                  {"build/sdpb", default_sdpb_args, args}, num_procs);
              }

            // Read c,B,y and check that (c - B.y) equals to the vector written to c_minus_By/c_minus_By.json
            check_c_minus_By(sdp_path, output_dir / "out", precision,
                             diff_precision,
                             runner.create_nested("check_c_minus_By"));

            // SDPB runs with --precision=<precision>
            // We check test output up to lower precision=<sdpb_output_diff_precision>
            // in order to neglect unimportant rounding errors
            diff_sdpb_output_dir(output_dir / "out", data_output_dir / "out",
                                 precision, diff_precision, sdpb_out_filenames,
                                 sdpb_out_txt_keys);
          }

          if(exists(data_output_dir / "spectrum.json"))
            {
              Named_Args_Map args{
                {"--pmpInfo", sdp_path + "/pmp_info.json"},
                {"--solution", (output_dir / "out").string()},
                {"--threshold", "1e-10"},
                {"--output", (output_dir / "spectrum.json").string()},
                {"--precision", std::to_string(precision)},
                {"--verbosity", "debug"},
              };
              runner.create_nested("spectrum")
                .mpi_run({"build/spectrum", args}, num_procs);

              // Cannot check block paths if the same spectrum.json
              // is generated several times by different PMP inputs,
              // e.g. pmp.xml and pmp.json.
              bool check_block_paths = pmp_paths.size() == 1;
              diff_spectrum(output_dir / "spectrum.json",
                            data_output_dir / "spectrum.json", precision,
                            diff_precision, check_block_paths);
            }
        }
    }
  };
}

namespace
{
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

  SECTION("1d")
  {
    INFO("SDPB test for a simple one-dimensional problem from SDPB Manual:");
    INFO("maximize (-y) s.t. (1 + x^4 + y * (x^4 / 12 + x^2)) >= 0)");

    End_To_End_Test test("1d");
    test.num_procs = 2;
    // Do not check normalization.json because it is absent in pmp-no-optional-fields.json
    test.check_sdp_normalization = false;
    // binary precision 664 is equivalent to decimal precision 200 used in SDPB.m
    // TODO generate input files with precision 768
    test.precision = 664;
    test.run();
  }

  SECTION("1d-old-sampling")
  {
    INFO("SDPB test for a simple one-dimensional problem from SDPB Manual:");
    INFO("maximize (-y) s.t. (1 + x^4 + y * (x^4 / 12 + x^2)) >= 0)");
    INFO("Use old sampling algorithm, sampling data specified explicitly in "
         "XML and JSON");

    End_To_End_Test test("1d-old-sampling");
    test.num_procs = 2;
    // Do not check normalization.json because it is absent for XML
    test.check_sdp_normalization = false;
    test.run();
  }

  SECTION("1d-preconditioning")
  {
    INFO("SDPB test for a simple one-dimensional problem from SDPB Manual:");
    INFO("maximize (-y) s.t. (1 + x^4 + y * (x^4 / 12 + x^2)) >= 0)");
    INFO("Same as 1d case, but the polynomial is multiplied by a custom "
         "preconditioning function f(x):");
    INFO("  pmp.json: f(x) = (0.2 + x)^0.3 * (10.1 + x + 3.1*x^2)^0.8");
    INFO("  pmp-const-10.json: f(x) = 10.0");
    INFO("Output is the same as in 1d, except for values in iterations.json "
         "and primal vector x_0.txt (its elements are divided by "
         "preconditioning values)");

    End_To_End_Test test("1d-preconditioning");

    test.num_procs = 2;
    test.precision = 664;
    // SDP is different because of preconditioning
    // TODO: compare everything except block_data?
    test.check_sdp = false;
    // Ignore iterations.json and x.txt since they depend on preconditioning
    test.sdpb_out_filenames = {"out.txt", "y.txt"};
    test.run();
  }

  SECTION("1d-duplicate-poles")
  {
    INFO("SDPB test for a simple one-dimensional problem from SDPB Manual:");
    INFO("maximize (-y) s.t. (1 + x^4 + y * (x^4 / 12 + x^2)) >= 0)");
    INFO("Here the matrix has a non-trivial prefactor with duplicate poles.");
    INFO("Primal objective, dual objective and the vector y should be similar "
         "to the simple 1d case. The vector x may differ.");
    End_To_End_Test test("1d-duplicate-poles");
    test.num_procs = 2;
    test.run();
  }

  SECTION("1d-constraints")
  {
    INFO("Same as 1d, but with nontrivial semidefiniteness constraints.");
    INFO("See mathematica/Tests.m");
    End_To_End_Test test("1d-constraints");
    test.num_procs = 2;
    test.run();
  }

  SECTION("1d-isolated-zeros")
  {
    INFO("maximize (-y) s.t. (1 + x^4 + y * (x^4 / 12 + x^2)) >= 0) for "
         "x=2/3, x=4/3, and x>=2");
    INFO("SDPB should find primal-dual optimal solution.");
    INFO("Spectrum should find isolated zero for the last block "
         "(corresponding to x=4/3).");
    End_To_End_Test test("1d-isolated-zeros");
    test.default_sdpb_args = boost::program_options::split_unix(
      "--checkpointInterval 3600 --maxRuntime 1340 "
      "--dualityGapThreshold 1.0e-30 --primalErrorThreshold 1.0e-30 "
      "--dualErrorThreshold 1.0e-30 --initialMatrixScalePrimal 1.0e20 "
      "--initialMatrixScaleDual 1.0e20 --feasibleCenteringParameter 0.1 "
      "--infeasibleCenteringParameter 0.3 --stepLengthReduction 0.7 "
      "--maxComplementarity 1.0e100 --maxIterations 1000 --verbosity 1 "
      "--procGranularity 1 --writeSolution x,y");
    test.num_procs = 1;
    test.check_sdp = false;
    // Write SDP to zip archive  to test that spectrum can read sdp.zip/pmp_info.json
    test.pmp2sdp_args = {{"--zip", ""}};
    test.run();
  }

  SECTION("dfibo-0-0-j=3-c=3.0000-d=3-s=6")
  {
    INFO("pmp2sdp+sdpb test for https://github.com/davidsd/sdpb/issues/124");
    INFO("sdp contains block with empty bilinear_bases_odd, "
         "which caused a bug.");
    INFO("Test data from Harvard cluster, gmp/6.2.1 mpfr/4.2.0");
    End_To_End_Test test("dfibo-0-0-j=3-c=3.0000-d=3-s=6");
    test.default_sdpb_args = boost::program_options::split_unix(
      "--findDualFeasible --findPrimalFeasible "
      "--initialMatrixScalePrimal 1e10 --initialMatrixScaleDual 1e10 "
      "--maxComplementarity 1e30 --dualErrorThreshold 1e-10 "
      "--primalErrorThreshold 1e-153 --maxRuntime 259200 "
      "--checkpointInterval 3600 --maxIterations 1000 "
      "--feasibleCenteringParameter=0.1 --infeasibleCenteringParameter=0.3 "
      "--stepLengthReduction=0.7 "
      "--maxSharedMemory=100K"); // forces split_factor=3 for Q window
    for(std::string sdp_format : {"", "bin", "json"})
      {
        DYNAMIC_SECTION(
          "format=" << (sdp_format.empty() ? "default(bin)" : sdp_format))
        {
          // write sdp to zip instead of plain directory
          bool zip = true;
          test.pmp2sdp_args = build_pmp2sdp_args(sdp_format, zip);
          test.run();
        }
      }
  }

  SECTION("SingletScalar_cT_test_nmax6")
  {
    INFO("SingletScalar_cT_test_nmax6 from "
         "https://gitlab.com/davidsd/scalars-3d/-/blob/master/src/Projects/"
         "Scalars3d/SingletScalar2020.hs");
    INFO("Test data is generated with SDPB 2.5.1 on Caltech cluster.");
    INFO("SDPB should find primal-dual optimal solution.");
    auto default_sdpb_args = boost::program_options::split_unix(
      "--checkpointInterval 3600 --maxRuntime 1340 "
      "--dualityGapThreshold 1.0e-30 --primalErrorThreshold 1.0e-30 "
      "--dualErrorThreshold 1.0e-30 --initialMatrixScalePrimal 1.0e20 "
      "--initialMatrixScaleDual 1.0e20 --feasibleCenteringParameter 0.1 "
      "--infeasibleCenteringParameter 0.3 --stepLengthReduction 0.7 "
      "--maxComplementarity 1.0e100 --maxIterations 1000 --verbosity 2 "
      "--procGranularity 1 --writeSolution x,y,z");
    SECTION("primal_dual_optimal")
    {
      End_To_End_Test test("SingletScalar_cT_test_nmax6/primal_dual_optimal");
      test.default_sdpb_args = default_sdpb_args;
      // This test is slow, we don't want to run it for both json and binary SDP.
      // json/bin correctness is checked by other tests,
      // so we use only binary SDP here

      test.check_sdp_normalization = true;
      // run_sdpb_twice=true to test checkpoint loading, see https://github.com/davidsd/sdpb/issues/219
      test.run_sdpb_twice = true;
      test.run();
    }

    SECTION("primal_dual_optimal_reduced")
    {
      INFO("Same as primal_dual_optimal, but with reducedPrefactor.");
      End_To_End_Test test(
        "SingletScalar_cT_test_nmax6/primal_dual_optimal_reduced");
      test.default_sdpb_args = default_sdpb_args;
      test.run();
    }
    SECTION("primal_dual_optimal_preconditioning")
    {
      INFO("Same as primal_dual_optimal, but with preconditioning.");
      End_To_End_Test test(
        "SingletScalar_cT_test_nmax6/primal_dual_optimal_preconditioning");
      test.default_sdpb_args = default_sdpb_args;
      test.check_sdp = false;
      test.run();
    }
  }

  SECTION("SingletScalarAllowed_test_nmax6")
  {
    INFO("SingletScalarAllowed_test_nmax6 from "
         "https://gitlab.com/davidsd/scalars-3d/-/blob/master/src/Projects/"
         "Scalars3d/SingletScalar2020.hs");
    INFO("Test data is generated with SDPB 2.5.1 on Caltech cluster.");
    std::string name = "SingletScalarAllowed_test_nmax6";
    auto default_sdpb_args = boost::program_options::split_unix(
      "--checkpointInterval 3600 --maxRuntime 1341 "
      "--dualityGapThreshold 1.0e-30 --primalErrorThreshold 1.0e-200 "
      "--dualErrorThreshold 1.0e-200 --initialMatrixScalePrimal 1.0e20 "
      "--initialMatrixScaleDual 1.0e20 --feasibleCenteringParameter 0.1 "
      "--infeasibleCenteringParameter 0.3 --stepLengthReduction 0.7 "
      "--maxComplementarity 1.0e100 --maxIterations 1000 --verbosity 2 "
      "--procGranularity 1 --writeSolution y,z "
      "--detectPrimalFeasibleJump --detectDualFeasibleJump "
      "--maxSharedMemory=100.1K"); // forces split_factor=3 for Q window; also test floating-point --maxSharedMemory value

    SECTION("primal_feasible_jump")
    {
      INFO("SDPB should detect primal feasible jump.");
      name += "/primal_feasible_jump";
      End_To_End_Test test(name);
      test.default_sdpb_args = default_sdpb_args;
      test.sdpb_out_txt_keys = {"terminateReason", "primalObjective",
                                "dualObjective", "dualityGap", "dualError"};
      test.run();
    }
    SECTION("dual_feasible_jump")
    {
      INFO("SDPB should detect dual feasible jump.");
      name += "/dual_feasible_jump";
      End_To_End_Test test(name);
      test.default_sdpb_args = default_sdpb_args;
      test.sdpb_out_txt_keys = {"terminateReason", "primalObjective",
                                "dualObjective", "dualityGap", "primalError"};

      test.run();
    }
  }
}
