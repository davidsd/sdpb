#include "integration_tests/common.hxx"

// GPU builds only, and it needs a GPU: sdpb --controlDir --gpu on the
// resident test's data. Its update_pmp and update_sdp change the spin-0 gap
// block, whose bilinear bases change with the gap, which --gpu doesn't take
// in resident mode: both must be answered with an error, the SDP in memory
// unchanged, and the session must go on (the solves around them succeed and
// quit ends it). Run on a GPU node with
//   ./build/integration_tests "resident_gpu" --mpirun="mpirun --oversubscribe"
#ifdef SDPB_GPU

#include <boost/program_options/parsers.hpp>

#include <fstream>
#include <sstream>

using namespace Test_Util;
using Named_Args_Map = Test_Case_Runner::Named_Args_Map;

namespace fs = std::filesystem;

namespace
{
  std::string read_text(const fs::path &path)
  {
    std::ifstream in(path, std::ios::binary);
    REQUIRE(in.good());
    std::stringstream ss;
    ss << in.rdbuf();
    return ss.str();
  }

  void write_text(const fs::path &path, const std::string &text)
  {
    std::ofstream out(path, std::ios::binary);
    out << text;
    REQUIRE(out.good());
  }
}

TEST_CASE("resident_gpu")
{
  INFO("sdpb --controlDir --gpu refuses a change of bilinear bases and goes on");
  const int precision = 768;
  const auto data_dir = Test_Config::test_data_dir / "end-to-end_tests"
                        / "SingletScalarAllowed_test_nmax6";
  const auto dir_A = data_dir / "dual_feasible_jump" / "input" / "json";
  const auto dir_B = data_dir / "primal_feasible_jump" / "input" / "json";
  const std::string changed_A
    = "mat_huHER-bgnrovCSCBgKZ40XuFTcLhpW6hNmHx0gQEbwM.json";
  const std::string changed_B
    = "mat_RfFLGRJJdlsbUjbyx3vPfOJ2OVffJDUPiaJFv6v7euI.json";

  Test_Case_Runner runner("resident_gpu");
  const auto &output_dir = runner.output_dir;
  fs::create_directories(output_dir);

  std::vector<std::string> common;
  for(const auto &entry : fs::directory_iterator(dir_A))
    if(entry.path().filename() != changed_A)
      common.push_back(fs::absolute(entry.path()).string());
  std::sort(common.begin(), common.end());
  const auto write_nsv = [&](const fs::path &path, const fs::path &first) {
    std::string text = fs::absolute(first).string() + '\0';
    for(const auto &file : common)
      text += file + '\0';
    write_text(path, text);
  };
  const auto nsv_A = output_dir / "pmp_A.nsv";
  const auto nsv_B = output_dir / "pmp_B.nsv";
  write_nsv(nsv_A, dir_A / changed_A);
  write_nsv(nsv_B, dir_B / changed_B);
  const auto sdp_A = output_dir / "sdp_A";
  const auto sdp_B = output_dir / "sdp_B";
  for(const auto &[nsv, sdp] :
      {std::pair{nsv_A, sdp_A}, std::pair{nsv_B, sdp_B}})
    {
      const Named_Args_Map args{{"--input", nsv.string()},
                                {"--output", sdp.string()},
                                {"--precision", std::to_string(precision)}};
      runner.create_nested("pmp2sdp_" + sdp.filename().string())
        .mpi_run({Test_Config::program("pmp2sdp"), args}, 1);
    }

  const auto sdpb_args = boost::program_options::split_unix(
    "--gpu --checkpointInterval 3600 --maxRuntime 1341 "
    "--dualityGapThreshold 1.0e-30 --primalErrorThreshold 1.0e-200 "
    "--dualErrorThreshold 1.0e-200 --initialMatrixScalePrimal 1.0e20 "
    "--initialMatrixScaleDual 1.0e20 --feasibleCenteringParameter 0.1 "
    "--infeasibleCenteringParameter 0.3 --stepLengthReduction 0.7 "
    "--maxComplementarity 1.0e100 --maxIterations 1000 --verbosity 1 "
    "--writeSolution y,z --detectPrimalFeasibleJump --detectDualFeasibleJump "
    "--precision "
    + std::to_string(precision));
  const auto ck = output_dir / "ck";
  const auto control = output_dir / "control";
  fs::create_directories(control);
  const auto solve = [&](int k) {
    return "{\"cmd\": \"solve\", \"args\": [\"--outDir="
           + (output_dir / ("out_" + std::to_string(k))).string() + "\"]}";
  };
  write_text(control / "cmd.1.json", solve(1));
  write_text(control / "cmd.2.json",
             "{\"cmd\": \"update_pmp\", \"pmp\": \"" + nsv_B.string()
               + "\"}");
  write_text(control / "cmd.3.json", solve(3));
  write_text(control / "cmd.4.json",
             "{\"cmd\": \"update_sdp\", \"sdp\": \"" + sdp_B.string()
               + "\"}");
  write_text(control / "cmd.5.json", solve(5));
  write_text(control / "cmd.6.json", "{\"cmd\": \"quit\"}");
  const Named_Args_Map resident_args{{"--sdpDir", sdp_A.string()},
                                     {"--checkpointDir", ck.string()},
                                     {"--controlDir", control.string()}};
  runner.create_nested("sdpb_resident")
    .mpi_run({Test_Config::program("sdpb"), sdpb_args, resident_args}, 1);
  REQUIRE(fs::exists(control / "exit.json"));
  for(int seq = 1; seq <= 6; ++seq)
    {
      CAPTURE(seq);
      const auto result
        = read_text(control / ("res." + std::to_string(seq) + ".json"));
      CAPTURE(result);
      if(seq == 2 || seq == 4)
        {
          REQUIRE(result.find("\"error\"") != std::string::npos);
          REQUIRE(result.find("bilinear bases") != std::string::npos);
        }
      else
        REQUIRE(result.find("\"error\"") == std::string::npos);
    }
  // The solves around the refused updates ran (on the SDP in memory, A)
  for(const int k : {1, 3, 5})
    REQUIRE(read_text(output_dir / ("out_" + std::to_string(k)) / "out.txt")
              .find("terminateReason") != std::string::npos);
}

#endif
