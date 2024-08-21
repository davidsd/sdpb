#include "integration_tests/common.hxx"
#include "integration_tests/util/Parse_Sdpb_Out_Txt.hxx"

#include <algorithm>
#include <map>
#include <boost/iostreams/filter/gzip.hpp>

using namespace Test_Util;
using namespace Test_Util::REQUIRE_Equal;
using Named_Args_Map = Test_Case_Runner::Named_Args_Map;

namespace fs = std::filesystem;

struct Sdplib_Case_Info
{
  std::string name;
  size_t m;
  size_t n;
  std::string optimal_binary_64;
  std::string optimal_gmp;

  void diff(const Float &obj) const
  {
    const auto optimal
      = optimal_gmp == "-" ? Float(optimal_binary_64) : Float(optimal_gmp);
  }
};

// Sort
struct Compare_Case_Group
{
  bool operator()(const std::vector<Sdplib_Case_Info> &lhs,
                  const std::vector<Sdplib_Case_Info> &rhs) const
  {
    return get_key(lhs) < get_key(rhs);
  }

private:
  static size_t get_key(const std::vector<Sdplib_Case_Info> &a)
  {
    size_t n_max = 0;
    for(const auto &item : a)
      {
        if(n_max < item.n)
          n_max = item.n;
      }
    return n_max;
  }
};

// These tests are slow and will be disabled by default, via "[.]"
TEST_CASE("SDPLIB", "[.]")
// TEST_CASE("SDPLIB")
{
  INFO("End-to-end tests for SDPA input files");
  INFO("Test cases from https://github.com/vsdp/SDPLIB");
  INFO("binary64 and high-precision (SDPA-GMP) objective values taken from "
       "https://github.com/nakatamaho/sdpa-gmp/blob/master/README.md");
  INFO("Note that sometimes binary64 and GMP results are very different, e.g. "
       "hinf1");
  int num_procs = 0;
  int precision = 200; //precision=200 fails with non-HPD on hinf7

  // Copied from https://github.com/nakatamaho/sdpa-gmp/blob/master/README.md
  std::vector<Sdplib_Case_Info> sdplib_cases{
    {"arch0", 174, 335, "5.66517e-01", "5.6651727321592959e-01"},
    {"arch2", 174, 335, "6.71515e-01", "6.7151540763990793e-01"},
    {"arch4", 174, 335, "9.726274e-01", "9.7262741740980893e-01"},
    {"arch8", 174, 335, "7.05698e+00", "7.0569800367002555e+00"},
    {"control1", 21, 15, "1.778463e+01", "1.7784626717523405e+01"},
    {"control2", 66, 30, "8.300000e+00", "8.2999999857902351e+00"},
    {"control3", 136, 45, "1.363327e+01", "1.3633266228377313e+01"},
    {"control4", 231, 60, "1.979423e+01", "1.9794230376537536e+01"},
    {"control5", 351, 75, "1.68836e+01", "1.6883599050955793e+01"},
    {"control6", 496, 90, "3.73044e+01", "3.7304412616333280e+01"},
    {"control7", 666, 105, "2.06251e+01", "2.0625072243801761e+01"},
    {"control8", 861, 120, "2.0286e+01", "2.0286360379531460e+01"},
    {"control9", 1081, 135, "1.46754e+01", "1.4675424692813939e+01"},
    {"control10", 1326, 150, "3.8533e+01", "3.8533032079581028e+01"},
    {"control11", 1596, 165, "3.1959e+01", "3.1958667450372498e+01"},
    {"eqaulG11", 801, 801, "6.291553e+02", "6.2915529282061428e+02"},
    {"equalG51", 1001, 1001, "4.005601e+03", "4.0056013154878550e+03"},
    {"gpp100", 101, 100, "-4.49435e+01", "-4.4943550775891146e+01"},
    {"gpp124-1", 125, 124, "-7.3431e+00", "-7.3430762652465377e+00"},
    {"gpp124-2", 125, 124, "-4.68623e+01", "-4.6862295072749908e+01"},
    {"gpp124-3", 125, 124, "-1.53014e+02", "-1.5301412730175306e+02"},
    {"gpp124-4", 125, 124, "-4.1899e+02", "-4.1898762587351130e+02"},
    {"gpp250-1", 250, 250, "-1.5445e+01", "-1.5444916882934067e+01"},
    {"gpp250-2", 250, 250, "-8.1869e+01", "-8.1868958840223643e+01"},
    {"gpp250-3", 250, 250, "-3.035e+02", "-3.0353932422884198e+02"},
    {"gpp250-4", 250, 250, "-7.473e+02", "-7.4732831101269269e+02"},
    {"gpp500-1", 501, 500, "-2.53e+01", "-2.5320543879075787e+01"},
    {"gpp500-2", 501, 500, "-1.5606e+02", "-1.5606038757941642e+02"},
    {"gpp500-3", 501, 500, "-5.1302e+02", "-5.1301760233182234e+02"},
    {"gpp500-4", 501, 500, "-1.56702e+03", "-1.5670187921561449e+03"},
    {"hinf1", 13, 14, "2.0326e+00", "2.6899806999311550e-05"},
    {"hinf2", 13, 16, "1.0967e+01", "1.0967055621049256e+01"},
    {"hinf3", 13, 16, "5.69e+01", "5.6940778009669388e+01"},
    {"hinf4", 13, 16, "2.74764e+02", "2.7149772617088246e+02"},
    {"hinf5", 13, 16, "3.63e+02", "3.3010372908768509e+02"},
    {"hinf6", 13, 16, "4.490e+02", "4.4892774532835125e+02"},
    {"hinf7", 13, 16, "3.91e+02", "1.5490470173390994e+02"},
    {"hinf8", 13, 16, "1.16e+02", "5.8449173799445524e+01"},
    {"hinf9", 13, 16, "2.3625e+02", "2.3624925825291886e+02"},
    {"hinf10", 21, 18, "1.09e+02", "1.3503382826378728e+01"},
    {"hinf11", 31, 22, "6.59e+01", "5.0979848897318626e+01"},
    {"hinf12", 43, 24, "2e-1", "1.5875842644916031e-13"},
    {"hinf13", 57, 30, "4.6e+01", "1.1816982963162776e-02"},
    {"hinf14", 73, 34, "1.30e+01", "2.8479790096215441e+00"},
    {"hinf15", 91, 37, "2.5e+01", "1.1224662755236866e-04"},
    {"infd1", 10, 30, "dual infeasible", "-1.1641795380064849e+05"},
    {"infd2", 10, 30, "dual infeasible", "-2.2856227327842922e+05"},
    {"infp1", 10, 30, "primal infeasible", "2.9943790805717583e+02"},
    {"infp2", 10, 30, "primal infeasible", "2.0749573188459872e+02"},
    {"maxG11", 800, 800, "6.291648e+02", "6.2916478300199902e+02"},
    {"maxG32", 2000, 2000, "1.567640e+03", "1.5676396446800114e+03"},
    {"maxG51", 1000, 1000, "4.003809e+03", "4.0062555216534127e+03"},
    {"maxG55", 5000, 5000, "9.999210e+03", "-"},
    {"maxG60", 7000, 7000, "1.522227e+04", "-"},
    {"mcp100", 100, 100, "2.261574e+02", "2.2615735148330884e+02"},
    {"mcp124-1", 124, 124, "1.419905e+02", "1.4199047709767370e+02"},
    {"mcp124-2", 124, 124, "2.698802e+02", "2.6988017064431990e+02"},
    {"mcp124-3", 124, 124, "4.677501e+02", "4.6775011428751099e+02"},
    {"mcp124-4", 124, 124, "8.644119e+02", "8.6441186405192719e+02"},
    {"mcp250-1", 250, 250, "3.172643e+02", "3.1726434034357982e+02"},
    {"mcp250-2", 250, 250, "5.319301e+02", "5.3193008393282009e+02"},
    {"mcp250-3", 250, 250, "9.811726e+02", "9.8117257166434770e+02"},
    {"mcp250-4", 250, 250, "1.681960e+03", "1.6819601121258921e+03"},
    {"mcp500-1", 500, 500, "5.981485e+02", "5.9814851691875962e+02"},
    {"mcp500-2", 500, 500, "1.070057e+03", "1.0700567662011862e+03"},
    {"mcp500-3", 500, 500, "1.847970e+03", "1.8479700215154574e+03"},
    {"mcp500-4", 500, 500, "3.566738e+03", "3.5667380499612209e+03"},
    {"qap5", 136, 26, "-4.360e+02", "-4.3600000000000000e+02"},
    {"qap6", 229, 37, "-3.8144e+02", "-3.8143840217367920e+02"},
    {"qap7", 358, 50, "-4.25e+02", "-4.2481971023200053e+02"},
    {"qap8", 529, 65, "-7.57e+02", "-7.5695525603861953e+02"},
    {"qap9", 748, 82, "-1.410e+03", "-1.4099410682401829e+03"},
    {"qap10", 1021, 101, "-1.093e+01", "-1.0926074684462390e+03"},
    {"qpG11", 800, 1600, "2.448659e+03", "2.4486591320079961e+03"},
    {"qpG51", 1000, 2000, "1.181000e+03", "1.1818000000000000e+04"},
    {"ss30", 132, 426, "2.02395e+01", "2.0239510569060567e+01"},
    {"theta1", 104, 50, "2.300000e+01", "2.3000000000000000e+01"},
    {"theta2", 498, 100, "3.287917e+01", "3.2879169015772581e+01"},
    {"theta3", 1106, 150, "4.216698e+01", "4.2166981488494406e+01"},
    {"theta4", 1949, 200, "5.032122e+01", "5.0321221951837344e+01"},
    {"theta5", 3028, 250, "5.723231e+01", "5.7232307282180003e+01"},
    {"theta6", 4375, 300, "6.347709e+01", "6.3477087177964743e+01"},
    {"thetaG11", 2401, 801, "4.000000e+02", "4.0000000000000000e+02"},
    {"thetaG51", 6910, 1001, "3.49000e+02", "3.4900000000000000e+02"},
    {"truss1", 6, 13, "-8.999996e+00", "-8.9999963152868905e+00"},
    {"truss2", 58, 133, "-1.233804e+02", "-1.2338035636407390e+02"},
    {"truss3", 27, 31, "-9.109996e+00", "-9.1099962092020534e+00"},
    {"truss4", 12, 19, "-9.009996e+00", "-9.0099962910045294e+00"},
    {"truss5", 208, 331, "-1.326357e+02", "-1.3263567797250604e+02"},
    {"truss6", 172, 451, "-9.01001e+02", "-9.0100140477088096e+02"},
    {"truss7", 86, 301, "-9.00001e+02", "-9.0000140369343463e+02"},
    {"truss8", 496, 628, "-1.331146e+02", "-1.3311458915226341e+02"}};

  // Group cases by common prefix (arch, control, equalG, ...)

  std::map<std::string, std::vector<Sdplib_Case_Info>> grouped_cases;
  for(const auto &sdplib_case : sdplib_cases)
    {
      const auto &name = sdplib_case.name;
      auto group = name;
      for(size_t i = 0; i < name.size(); i++)
        {
          if(!isalpha(name[i]))
            {
              group = name.substr(0, i);
              break;
            }
        }
      grouped_cases[group].push_back(sdplib_case);
    }

  // Sort groups, so that we start with the fastest tests

  std::vector<std::string> sorted_group_keys;
  for(const auto &[k, _] : grouped_cases)
    {
      sorted_group_keys.push_back(k);
    }

  // Take the largest size estimate for all test cases.
  // The size is estimated as the total size of all SDP matrices F_0, F_1,..F_m
  // There are (m + 1) matrices, each has the same size (n x n).
  const auto get_size_estimate = [](const std::vector<Sdplib_Case_Info> &a) {
    size_t size = 0;
    for(const auto &item : a)
      {
        const auto primal_dim = item.m;
        const auto sdp_dim = item.n;
        size = std::max(size, (primal_dim + 1) * sdp_dim * sdp_dim);
      }
    return size;
  };

  const auto compare_groups
    = [&get_size_estimate](const std::vector<Sdplib_Case_Info> &lhs,
                           const std::vector<Sdplib_Case_Info> &rhs) {
        return get_size_estimate(lhs) < get_size_estimate(rhs);
      };
  const auto compare_group_keys
    = [&compare_groups, &grouped_cases](const std::string &lhs,
                                        const std::string &rhs) {
        return compare_groups(grouped_cases.at(lhs), grouped_cases.at(rhs));
      };

  std::sort(sorted_group_keys.begin(), sorted_group_keys.end(),
            compare_group_keys);

  // Run tests
  for(const auto &group_key : sorted_group_keys)
    {
      DYNAMIC_SECTION(group_key)
      {
        for(auto &[name, primal_dim_m, sdp_dim_n, bin64_objective,
                   sdpa_gmp_objective] : grouped_cases.at(group_key))
          {
            DYNAMIC_SECTION(name)
            {
              CAPTURE(name);
              CAPTURE(primal_dim_m);
              CAPTURE(sdp_dim_n);
              CAPTURE(bin64_objective);
              CAPTURE(sdpa_gmp_objective);

              Test_Case_Runner runner("SDPLIB/" + name);
              const auto &output_dir = runner.output_dir;

              //sdpb
              {
                const auto sdp_path = Test_Config::test_data_dir
                                      / "SDPLIB/data" / (name + ".dat-s");
                // Default parameters used in sdpa-gmp
                std::string default_sdpb_args
                  = "--dualityGapThreshold 1.0e-30 "
                    "--initialMatrixScalePrimal 1.0e+4 "
                    "--initialMatrixScaleDual 1.0e+4 "
                    "--maxComplementarity 4.0e+8 "
                    "--feasibleCenteringParameter 0.3 stepLengthReduction 0.9 "
                    "--primalErrorThreshold 1.0e-30 "
                    "--dualErrorThreshold 1.0e-30 "
                    "--writeSolution x --noFinalCheckpoint";

                Named_Args_Map args{
                  {"--precision", std::to_string(precision)},
                  {"--sdpDir", sdp_path},
                  {"--outDir", (output_dir / "out").string()},
                  {"--checkpointDir", (output_dir / "ck").string()}};

                runner.create_nested("sdpb").mpi_run(
                  {"build/sdpb", default_sdpb_args}, args, num_procs);
              }
              //Check output
              {
                Float_Binary_Precision fbp(precision, precision);

                const auto out_txt_path = output_dir / "out/out.txt";
                INFO(out_txt_path.string());
                Parse_Sdpb_Out_Txt out_txt(out_txt_path);

                CAPTURE(out_txt.terminate_reason);
                // TODO some tests are infeasible
                REQUIRE(out_txt.terminate_reason
                        == "found primal-dual optimal solution");
                const auto &primal_objective
                  = out_txt.float_map["primalObjective"];
                CAPTURE(primal_objective);

                const auto &expected_objective_str = sdpa_gmp_objective == "-"
                                                       ? bin64_objective
                                                       : sdpa_gmp_objective;
                const Float expected_objective(expected_objective_str);

                // How many digits for expected_objective are known
                size_t digits10 = 0;
                for(auto c : expected_objective_str)
                  {
                    if(std::isdigit(c))
                      ++digits10;
                    // stop at exponent
                    if(c == 'e')
                      break;
                  }

                CAPTURE(digits10);
                // convert to binary digits
                const int diff_prec
                  = std::floor((digits10 - 1) * std::log2(10.0));
                DIFF_PREC(primal_objective, expected_objective, diff_prec);
              }
            }
          }
      }
    }
}
