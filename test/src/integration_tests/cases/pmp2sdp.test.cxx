#include "integration_tests/common.hxx"

#include <algorithm>

#include <filesystem>
#include <fstream>
#include <iterator>
#include <map>

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

  SECTION("number of ranks")
  {
    INFO("pmp2sdp on 1 and 3 ranks writes the same SDP, byte for byte "
         "(except control.json, which records the command). pmp_info.json "
         "used to depend on which rank owned each block: its prefactors "
         "were sent through El::BigFloat, which rounds Boost_Float values "
         "at 768 bits.");
    const auto input = Test_Config::test_data_dir / "end-to-end_tests"
                       / "SingletScalar_cT_test_nmax6"
                       / "primal_dual_optimal" / "input" / "pmp.nsv";
    Test_Util::Test_Case_Runner runner("pmp2sdp/number_of_ranks");
    Test_Util::Test_Case_Runner::Named_Args_Map args(default_args);
    args["--input"] = input.string();
    args["--precision"] = "768";
    std::map<int, fs::path> sdp;
    for(const int num_procs : {1, 3})
      {
        sdp[num_procs] = runner.output_dir / ("sdp_np" + std::to_string(num_procs));
        args["--output"] = sdp[num_procs].string();
        runner.create_nested("np" + std::to_string(num_procs))
          .mpi_run({Test_Config::program("pmp2sdp"), args}, num_procs);
      }
    auto read = [](const fs::path &path) {
      std::ifstream is(path, std::ios::binary);
      REQUIRE(is.good());
      return std::string(std::istreambuf_iterator<char>(is), {});
    };
    size_t num_files = 0;
    for(const auto &entry : fs::directory_iterator(sdp[1]))
      {
        const auto name = entry.path().filename();
        if(name == "control.json")
          continue;
        CAPTURE(name);
        REQUIRE(fs::exists(sdp[3] / name));
        REQUIRE(read(entry.path()) == read(sdp[3] / name));
        ++num_files;
      }
    REQUIRE(num_files > 2);
  }

  SECTION("reuseSdp")
  {
    INFO("pmp2sdp --reuseSdp: in a list of input files where some have new "
         "paths (as in a scan, where only some files change), the blocks of "
         "the others are linked from the previous SDP. The result must be "
         "byte-identical to a full build of the same list, except "
         "control.json, which records the command.");
    const auto orig_input
      = Test_Config::test_data_dir / "end-to-end_tests"
        / "SingletScalar_cT_test_nmax6" / "primal_dual_optimal" / "input";
    Test_Util::Test_Case_Runner runner("pmp2sdp/reuseSdp");
    const auto input = runner.output_dir / "input";
    fs::create_directories(runner.output_dir);
    fs::copy(orig_input, input, fs::copy_options::recursive);

    // the list of input files (NUL-separated, relative to the list)
    std::vector<std::string> files;
    {
      std::ifstream is(input / "pmp.nsv", std::ios::binary);
      std::string entry;
      while(std::getline(is, entry, '\0'))
        if(!entry.empty())
          files.push_back(entry);
    }
    REQUIRE(files.size() > 3);
    // the objective, the normalization and one more matrix get new paths
    // (same contents), as a changed file would in a scan; the list is also
    // reversed, so that the reused blocks get new block indices
    std::vector<std::string> new_files(files);
    for(size_t i = 0; i < 3; ++i)
      {
        new_files[i] = "json/new_" + std::to_string(i) + ".json";
        fs::copy_file(input / files[i], input / new_files[i]);
      }
    std::reverse(new_files.begin(), new_files.end());
    {
      std::ofstream os(input / "pmp_new.nsv", std::ios::binary);
      for(const auto &f : new_files)
        os << f << '\0';
    }

    Test_Util::Test_Case_Runner::Named_Args_Map args(default_args);
    args["--precision"] = "768";
    args["--input"] = (input / "pmp.nsv").string();
    const auto sdp_old = runner.output_dir / "sdp_old";
    args["--output"] = sdp_old.string();
    runner.create_nested("old").mpi_run({Test_Config::program("pmp2sdp"), args});

    args["--input"] = (input / "pmp_new.nsv").string();
    const auto sdp_full = runner.output_dir / "sdp_full";
    args["--output"] = sdp_full.string();
    runner.create_nested("full").mpi_run({Test_Config::program("pmp2sdp"), args});

    Test_Util::Test_Case_Runner::Named_Args_Map reuse_args(args);
    reuse_args["--reuseSdp"] = sdp_old.string();
    const auto sdp_reuse = runner.output_dir / "sdp_reuse";
    reuse_args["--output"] = sdp_reuse.string();
    runner.create_nested("reuse").mpi_run({Test_Config::program("pmp2sdp"), reuse_args});

    auto read = [](const fs::path &path) {
      std::ifstream is(path, std::ios::binary);
      return std::string(std::istreambuf_iterator<char>(is), {});
    };
    size_t num_files = 0, num_linked = 0;
    for(const auto &entry : fs::directory_iterator(sdp_full))
      {
        const auto name = entry.path().filename();
        if(name == "control.json")
          continue;
        CAPTURE(name);
        REQUIRE(fs::exists(sdp_reuse / name));
        REQUIRE(read(sdp_reuse / name) == read(entry.path()));
        ++num_files;
        if(fs::hard_link_count(sdp_reuse / name) > 1)
          ++num_linked;
      }
    for(const auto &entry : fs::directory_iterator(sdp_reuse))
      REQUIRE(fs::exists(sdp_full / entry.path().filename()));
    REQUIRE(num_files > 0);
    // the unchanged files' blocks are linked, not rebuilt
    REQUIRE(num_linked > 0);

    // a previous SDP at another precision must be refused
    {
      auto refused_args = reuse_args;
      refused_args["--precision"] = "640";
      refused_args["--output"] = (runner.output_dir / "sdp_refused").string();
      runner.create_nested("refused").mpi_run(
        {Test_Config::program("pmp2sdp"), refused_args}, 2, 1,
        "requires the old SDP to have the same precision");
    }
    // ... and so must a changed normalization (read from a new file)
    {
      const auto normalization_file = std::find_if(
        new_files.begin(), new_files.end(), [&](const std::string &f) {
          return f.find("new_") != std::string::npos
                 && read(input / f).find("normalization")
                      != std::string::npos;
        });
      REQUIRE(normalization_file != new_files.end());
      std::string text = read(input / *normalization_file);
      // change the first digit after the first '"' of a number
      const auto pos = text.find_first_of("123456789",
                                          text.find("normalization"));
      REQUIRE(pos != std::string::npos);
      text[pos] = text[pos] == '9' ? '8' : char(text[pos] + 1);
      {
        std::ofstream os(input / *normalization_file, std::ios::binary);
        os << text;
      }
      auto refused_args = reuse_args;
      refused_args["--output"]
        = (runner.output_dir / "sdp_refused_normalization").string();
      runner.create_nested("refused_normalization")
        .mpi_run({Test_Config::program("pmp2sdp"), refused_args}, 2, 1,
                 "the normalization differs");
    }
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
