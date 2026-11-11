#include "integration_tests/common.hxx"

#include <fstream>
#include <map>
#include <rapidjson/document.h>

#include <filesystem>

namespace fs = std::filesystem;
using namespace std::string_literals;

TEST_CASE("pmp2sdp")
{
  unsigned int precision = 512;
  unsigned int diff_precision = 392;
  Test_Util::Test_Case_Runner::Named_Args_Map default_args{
    {"--precision", std::to_string(precision)},
    {"--verbosity", "debug"},
  };

  SECTION("run")
  {
    std::string input_format = GENERATE("json", "m");
    DYNAMIC_SECTION(input_format)
    {
      auto data_dir = Test_Config::test_data_dir / "pmp2sdp" / input_format;
      auto sdp_orig = data_dir / "sdp_orig";

      for(const std::string &input_filename :
          {"pmp."s + input_format, "file_list.nsv"s})
        DYNAMIC_SECTION(input_filename)
        {
          Test_Util::Test_Case_Runner runner("pmp2sdp/"s + input_format + "/"
                                             + input_filename);
          Test_Util::Test_Case_Runner::Named_Args_Map args(default_args);
          args["--input"] = (data_dir / input_filename).string();
          auto sdp_path = (runner.output_dir / "sdp").string();
          args["--output"] = sdp_path;

          runner.create_nested("run").mpi_run(
            {Test_Config::program("pmp2sdp"), args});

          Test_Util::REQUIRE_Equal::diff_sdp(sdp_path, sdp_orig, precision,
                                             diff_precision,
                                             runner.create_nested("diff"));

          REQUIRE(fs::file_size(sdp_path + ".profiling/profiling.0") > 0);
          REQUIRE(fs::file_size(sdp_path + ".profiling/profiling.1") > 0);
        }
    }
  }

  SECTION("outputFormat")
  {
    INFO("Check different --outputFormat options");
    // TODO use input from "1d" and update output
    auto data_dir
      = Test_Config::test_data_dir / "end-to-end_tests" / "1d-old-sampling";

    for(std::string output_format : {"", "bin", "json"})
      {
        auto format_description
          = output_format.empty() ? "default(bin)" : output_format;
        CAPTURE(output_format);
        CAPTURE(format_description);

        DYNAMIC_SECTION(format_description)
        {
          Test_Util::Test_Case_Runner runner("pmp2sdp/outputFormat"s + "/"
                                             + format_description);
          Test_Util::Test_Case_Runner::Named_Args_Map args(default_args);
          args["--input"] = (data_dir / "input" / "pmp.json").string();
          auto sdp_path = runner.output_dir / "sdp";
          args["--output"] = sdp_path.string();
          if(!output_format.empty())
            args["--outputFormat"] = output_format;
          runner.create_nested("run").mpi_run(
            {Test_Config::program("pmp2sdp"), args});

          {
            INFO("Check that pmp2sdp actually uses --outputFormat="
                 << format_description);
            auto block_data_0_path
              = sdp_path
                / ("block_data_0."
                   + (output_format.empty() ? "bin" : output_format));
            CAPTURE(block_data_0_path);
            REQUIRE(is_regular_file(block_data_0_path));
          }

          auto sdp_orig = data_dir / "output" / "sdp";

          Test_Util::REQUIRE_Equal::diff_sdp(sdp_path, sdp_orig, precision,
                                             diff_precision,
                                             runner.create_nested("diff"));
        }
      }
  }

  SECTION("reducedPrefactorNodesAtInfinity")
  {
    INFO("Every polynomial is P = Pi_dropped * S with deg S = d, so the "
         "reduced numerator P / Pi_dropped is a polynomial of degree d and "
         "conditions at infinity must give the same SDP as plain "
         "interpolation.");
    auto data_dir = Test_Config::test_data_dir / "pmp2sdp" / "nodes-at-infinity";
    auto input = (data_dir / "pmp.json").string();

    Test_Util::Test_Case_Runner runner("pmp2sdp/nodes-at-infinity");
    Test_Util::Test_Case_Runner::Named_Args_Map args(default_args);
    args["--input"] = input;
    args["--maxNumPoles"] = "2";
    args["--outputFormat"] = "json";
    const auto plain_sdp = runner.output_dir / "sdp-plain";
    args["--output"] = plain_sdp.string();
    runner.create_nested("plain").mpi_run({Test_Config::program("pmp2sdp"), args});

    for(const std::string num_nodes : {"1", "2", "5"})
      DYNAMIC_SECTION("M=" << num_nodes)
      {
        auto nested = runner.create_nested("M=" + num_nodes);
        Test_Util::Test_Case_Runner::Named_Args_Map args_inf(args);
        args_inf["--reducedPrefactorNodesAtInfinity"] = num_nodes;
        const auto sdp = nested.output_dir / "sdp";
        args_inf["--output"] = sdp.string();
        nested.create_nested("run").mpi_run({Test_Config::program("pmp2sdp"), args_inf});
        Test_Util::REQUIRE_Equal::diff_sdp(sdp, plain_sdp, precision,
                                           diff_precision,
                                           nested.create_nested("diff"));
      }
  }

  SECTION("reducedPrefactorNodesAtInfinity on a generic PMP")
  {
    INFO("On SingletScalar_cT_test_nmax6 with --maxNumPoles 14 (blocks have "
         "11 to 17 poles), the nodes at infinity must change the blocks "
         "whose poles were reduced, and leave the others byte-identical.");
    const auto input = Test_Config::test_data_dir / "end-to-end_tests"
                       / "SingletScalar_cT_test_nmax6"
                       / "primal_dual_optimal" / "input" / "pmp.nsv";
    Test_Util::Test_Case_Runner runner("pmp2sdp/nodes-at-infinity-generic");
    Test_Util::Test_Case_Runner::Named_Args_Map args(default_args);
    args["--input"] = input.string();
    args["--precision"] = "768";
    args["--maxNumPoles"] = "14";
    const auto plain_sdp = runner.output_dir / "sdp-plain";
    args["--output"] = plain_sdp.string();
    runner.create_nested("plain").mpi_run({Test_Config::program("pmp2sdp"), args});
    args["--reducedPrefactorNodesAtInfinity"] = "2";
    const auto inf_sdp = runner.output_dir / "sdp-inf";
    args["--output"] = inf_sdp.string();
    runner.create_nested("M=2").mpi_run({Test_Config::program("pmp2sdp"), args});

    auto read = [](const fs::path &path) {
      std::ifstream is(path, std::ios::binary);
      REQUIRE(is.good());
      return std::string(std::istreambuf_iterator<char>(is), {});
    };
    // A block changes iff its poles were reduced (pmp_info.json: fewer
    // reducedPrefactor poles than prefactor poles)
    std::map<size_t, bool> reduced;
    {
      rapidjson::Document doc;
      const auto text = read(plain_sdp / "pmp_info.json");
      doc.Parse(text.c_str());
      REQUIRE(!doc.HasParseError());
      REQUIRE(doc.IsArray());
      for(const auto &block : doc.GetArray())
        reduced[block["index"].GetUint64()]
          = block["reducedPrefactor"]["poles"].Size()
            < block["prefactor"]["poles"].Size();
    }
    size_t num_changed = 0, num_same = 0, num_reduced = 0;
    for(const auto &[index, is_reduced] : reduced)
      num_reduced += is_reduced;
    for(const auto &entry : fs::directory_iterator(plain_sdp))
      {
        const auto name = entry.path().filename().string();
        CAPTURE(name);
        if(name == "control.json")
          continue;
        const bool same = read(entry.path()) == read(inf_sdp / name);
        if(name.rfind("block_data_", 0) == 0)
          {
            const size_t index = std::stoul(name.substr(11));
            REQUIRE(reduced.count(index) == 1);
            REQUIRE(same == !reduced.at(index));
            (same ? num_same : num_changed) += 1;
          }
        else if(name != "pmp_info.json")
          // block_info_*, objectives, normalization: unchanged
          REQUIRE(same);
      }
    CAPTURE(num_changed);
    CAPTURE(num_same);
    REQUIRE(num_changed == num_reduced);
    REQUIRE(num_changed > 0);
    REQUIRE(num_same > 0);
  }

  SECTION("filesystem errors")
  {
    INFO("pmp2sdp should fail due to invalid input/output arguments");
    auto data_dir = Test_Config::test_data_dir / "end-to-end_tests" / "1d";
    int num_procs = 2;

    auto input = (data_dir / "input" / "pmp.json").string();

    SECTION("sdp_temp_not_empty")
    {
      const Test_Util::Test_Case_Runner runner("pmp2sdp/sdp_temp_not_empty");
      fs::create_directories(runner.output_dir);

      auto sdp_dir = runner.output_dir / "sdp";

      {
        auto sdp_temp_dir = runner.output_dir / "sdp_temp";
        fs::create_directories(runner.output_dir / "sdp_temp");
        std::ofstream os(sdp_temp_dir / "some_file.json");
        os << "some data";
      }

      Test_Util::Test_Case_Runner::Named_Args_Map args(default_args);
      args["--input"] = input;
      args["--output"] = sdp_dir.string();
      // We allow pmp2sdp to overwrite exsiting sdp_temp
      runner.mpi_run({Test_Config::program("pmp2sdp"), args}, num_procs, 0,
                     "exists and will be overwritten");
    }
    SECTION("sdp_dir_exists")
    {
      INFO("Prohibit writing to existing directory:");
      const Test_Util::Test_Case_Runner runner("pmp2sdp/sdp_dir_exists");
      fs::create_directories(runner.output_dir);

      auto sdp_dir = runner.output_dir / "sdp";
      {
        fs::create_directories(sdp_dir);
        std::ofstream os(sdp_dir / "control.json");
        os << "some data";
      }

      Test_Util::Test_Case_Runner::Named_Args_Map args(default_args);
      args["--input"] = input;
      args["--output"] = sdp_dir.string();
      // We allow pmp2sdp to overwrite exsiting sdp
      runner.mpi_run({Test_Config::program("pmp2sdp"), args}, num_procs, 0,
                     "exists and will be overwritten");
    }

    SECTION("cannot_write_zip")
    {
      INFO("Prohibit output sdp.zip from writing and check that pmp2sdp "
           "fails:");
      const Test_Util::Test_Case_Runner runner("pmp2sdp/cannot_write_zip");
      fs::create_directories(runner.output_dir);

      auto sdp_readonly_zip = runner.output_dir / "sdp.readonly.zip";
      std::ofstream os(sdp_readonly_zip);
      os << "INVALID ZIP";
      permissions(sdp_readonly_zip, fs::perms::others_read);

      Test_Util::Test_Case_Runner::Named_Args_Map args(default_args);
      args["--input"] = input;
      args["--output"] = sdp_readonly_zip.string();
      args["--zip"] = "";
      // We allow pmp2sdp to overwrite exsiting sdp
      runner.mpi_run({Test_Config::program("pmp2sdp"), args}, num_procs, 0,
                     "exists and will be overwritten");
    }

    SECTION("invalid_nsv")
    {
      INFO("pmp2sdp should fail on an incorrect .nsv");
      const Test_Util::Test_Case_Runner runner("pmp2sdp/invalid_nsv");
      fs::create_directories(runner.output_dir);
      auto invalid_nsv
        = (runner.output_dir / "invalid_file_list.nsv").string();
      {
        std::ofstream os(invalid_nsv);
        os << "no_such_file.xml";
      }

      Test_Util::Test_Case_Runner::Named_Args_Map args(default_args);
      args["--input"] = invalid_nsv;
      auto sdp_path = (runner.output_dir / "sdp").string();
      args["--output"] = sdp_path;

      auto sdp_invalid_zip = (runner.output_dir / "sdp.invalid.zip").string();
      runner.mpi_run({Test_Config::program("pmp2sdp"), args}, num_procs, 1,
                     "No such file");
    }
  }

  SECTION("duplicate objectives")
  {
    INFO("Two PMP files having same objective vector.");
    INFO("pmp2sdp should not fail, see "
         "https://github.com/davidsd/sdpb/issues/251");
    INFO("Same input as end-to-end_tests/1d-constraints/pmp.xml, but split "
         "into two files");

    auto data_dir
      = Test_Config::test_data_dir / "pmp2sdp" / "duplicate-objectives";
    int num_procs = GENERATE(1, 2);
    auto input_nsv = (data_dir / "pmp.nsv").string();

    DYNAMIC_SECTION("num_procs=" << num_procs)
    {
      const Test_Util::Test_Case_Runner runner(
        "pmp2sdp/duplicate-objectives/num_procs=" + std::to_string(num_procs));

      Test_Util::Test_Case_Runner::Named_Args_Map args(default_args);
      args["--input"] = input_nsv;
      auto sdp_path = (runner.output_dir / "sdp").string();
      args["--output"] = sdp_path;

      runner.create_nested("run").mpi_run(
        {Test_Config::program("pmp2sdp"), args}, num_procs);

      auto sdp_orig = Test_Config::test_data_dir / "end-to-end_tests"
                      / "1d-constraints" / "output" / "sdp";
      Test_Util::REQUIRE_Equal::diff_sdp(sdp_path, sdp_orig, precision,
                                         diff_precision,
                                         runner.create_nested("diff"));
    }
  }

  SECTION("objectives are different")
  {
    INFO("pmp2sdp should fail because two files have different objectives");

    auto data_dir
      = Test_Config::test_data_dir / "pmp2sdp" / "wrong-objectives";
    int num_procs = GENERATE(1, 2);
    auto input_nsv = (data_dir / "pmp.nsv").string();

    DYNAMIC_SECTION("num_procs=" << num_procs)
    {
      const Test_Util::Test_Case_Runner runner(
        "pmp2sdp/wrong-objectives/num_procs=" + std::to_string(num_procs));

      Test_Util::Test_Case_Runner::Named_Args_Map args(default_args);
      args["--input"] = input_nsv;
      auto sdp_path = (runner.output_dir / "sdp").string();
      args["--output"] = sdp_path;
      runner.mpi_run({Test_Config::program("pmp2sdp"), args}, num_procs, 1,
                     "Found different objective vectors in input files:");
    }
  }
}
