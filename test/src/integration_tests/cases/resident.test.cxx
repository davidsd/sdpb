#include "integration_tests/common.hxx"

#include <boost/program_options/parsers.hpp>

#include <rapidjson/document.h>

#include <chrono>
#include <condition_variable>
#include <fstream>
#include <iostream>
#include <map>
#include <mutex>
#include <set>
#include <sstream>
#include <thread>

// Resident mode (sdpb --controlDir, the same sequence of SDPs
// solved by one resident process gives the same results, bit for bit, as
// sdpb runs chained through a checkpoint directory.
//
// Data: SingletScalarAllowed_test_nmax6. Its dual_feasible_jump and
// primal_feasible_jump inputs differ in one block (the spin-0 gap block), so
// solving A, B, A switches between excluded and allowed. The changed file
// comes first (block 0), and block_timings makes block 0 cost twice all the
// others together, so that with 4 ranks it is owned by a group of 2 ranks
// (the other blocks by two single-rank groups): replacing it in memory then
// goes through copy_matrix_from_root on a multi-rank grid.
//
// File route: pmp2sdp (A and B), then sdpb per solve with one checkpoint dir.
// Resident: one sdpb --controlDir that reads pre-written commands:
//   solve, save_iterate (in memory and as a checkpoint); update_pmp B, solve;
//   update_sdp A, solve; then update_sdp B and solve B three more times, from
//   the iterate saved after the first solve: restore_iterate by name,
//   restore_iterate from the saved checkpoint, and "start": {"iterate"};
//   quit.
// Compared after each solve: out.txt (except the runtime), y.txt, z.txt and
// every checkpoint file. The three restored solves stop after
// restore_iterations iterations, and must each equal the file route's solve
// of B from the checkpoint of the first solve, stopped there too: a restore
// that misses any of x, X, y, Y starts B from part of the third solve's
// iterate, and differs from the first step on.
//
// Both tests print a line per step, with the elapsed time, and one every 30 s
// while a step runs: CI services stop a job that is silent for minutes.

using namespace Test_Util;
using Named_Args_Map = Test_Case_Runner::Named_Args_Map;

namespace fs = std::filesystem;

namespace
{
  std::string read_file(const fs::path &path)
  {
    std::ifstream in(path, std::ios::binary);
    REQUIRE(in.good());
    std::stringstream ss;
    ss << in.rdbuf();
    return ss.str();
  }

  void write_file(const fs::path &path, const std::string &text)
  {
    std::ofstream out(path, std::ios::binary);
    out << text;
    REQUIRE(out.good());
  }

  std::string out_txt_without_runtime(const fs::path &out_dir)
  {
    std::istringstream in(read_file(out_dir / "out.txt"));
    std::string line, result;
    while(std::getline(in, line))
      if(line.rfind("Solver runtime", 0) != 0)
        result += line + "\n";
    return result;
  }

  std::string terminate_reason(const fs::path &out_dir)
  {
    std::istringstream in(read_file(out_dir / "out.txt"));
    std::string line;
    std::getline(in, line);
    return line;
  }

  // Copy the newest checkpoint generation of ck to dest
  void snapshot_checkpoint(const fs::path &ck, const fs::path &dest)
  {
    fs::create_directories(dest);
    for(const auto &entry : fs::directory_iterator(ck))
      if(entry.path().filename().string().rfind("checkpoint_", 0) == 0)
        fs::copy_file(entry.path(), dest / entry.path().filename(),
                      fs::copy_options::overwrite_existing);
  }

  // The files of checkpoint generation gen in dir (with gen < 0: of the only
  // generation there), by the name after "checkpoint_<generation>_", one per
  // rank
  std::map<std::string, fs::path>
  checkpoint_files_by_rank(const fs::path &dir, const int gen = -1)
  {
    std::map<std::string, fs::path> result;
    std::set<std::string> generations;
    for(const auto &entry : fs::directory_iterator(dir))
      {
        const auto name = entry.path().filename().string();
        if(name.rfind("checkpoint_", 0) != 0)
          continue;
        const auto rest = name.substr(std::string("checkpoint_").size());
        const auto pos = rest.find('_');
        REQUIRE(pos != std::string::npos);
        if(gen >= 0 && rest.substr(0, pos) != std::to_string(gen))
          continue;
        generations.insert(rest.substr(0, pos));
        result[rest.substr(pos + 1)] = entry.path();
      }
    REQUIRE(generations.size() == 1);
    return result;
  }

  std::vector<std::string> checkpoint_files(const fs::path &dir, int gen)
  {
    std::vector<std::string> result;
    for(const auto &entry : fs::directory_iterator(dir))
      {
        const auto name = entry.path().filename().string();
        if(name.rfind("checkpoint_" + std::to_string(gen) + "_", 0) == 0)
          result.push_back(name);
      }
    std::sort(result.begin(), result.end());
    return result;
  }

  // Block 0 costs twice all the other blocks together: with 4 ranks,
  // compute_block_grid_mapping gives it floor(4 * 2/3) = 2 ranks
  std::string block_timings_for(const size_t num_blocks)
  {
    std::string result = std::to_string(2 * (num_blocks - 1)) + "\n";
    for(size_t i = 1; i < num_blocks; ++i)
      result += "1\n";
    return result;
  }

  // Prints each step with the time since the test started, and the current
  // step again every 30 s while it runs. While a control directory is
  // watched, also prints each result the resident session writes there.
  class Progress
  {
  public:
    explicit Progress(std::string name) : name(std::move(name)) {}
    ~Progress()
    {
      {
        std::lock_guard lock(mutex);
        stop = true;
      }
      cv.notify_all();
      watcher.join();
    }
    Progress(const Progress &) = delete;
    Progress &operator=(const Progress &) = delete;

    void step(const std::string &what)
    {
      std::lock_guard lock(mutex);
      print_results();
      current = what;
      print(what);
    }
    // Report res.1.json, res.2.json, ... in dir as they appear
    void watch(const fs::path &dir)
    {
      std::lock_guard lock(mutex);
      control_dir = dir;
      next_result = 1;
    }

  private:
    using Clock = std::chrono::steady_clock;
    const std::string name;
    const Clock::time_point start = Clock::now();
    Clock::time_point last_print = start;
    std::mutex mutex;
    std::condition_variable cv;
    bool stop = false;
    std::string current;
    fs::path control_dir;
    int next_result = 1;
    // Last: it starts once everything above is initialized
    std::thread watcher{[this] { run(); }};

    // With the mutex held
    void print(const std::string &what)
    {
      last_print = Clock::now();
      std::cout << name << " ["
                << std::chrono::duration_cast<std::chrono::seconds>(last_print
                                                                    - start)
                     .count()
                << " s] " << what << std::endl;
    }
    // With the mutex held
    void print_results()
    {
      std::error_code ec;
      while(!control_dir.empty()
            && fs::exists(control_dir
                            / ("res." + std::to_string(next_result) + ".json"),
                          ec))
        print("resident session: answered cmd." + std::to_string(next_result++)
              + ".json");
    }
    void run()
    {
      std::unique_lock lock(mutex);
      while(
        !cv.wait_for(lock, std::chrono::seconds(1), [this] { return stop; }))
        {
          print_results();
          if(Clock::now() - last_print >= std::chrono::seconds(30))
            print("still running: " + current);
        }
    }
  };
}

TEST_CASE("resident")
{
  INFO("sdpb --controlDir reproduces chained sdpb runs bitwise");
  const int num_procs = 4;
  const int precision = 384;
  // Iterations of each solve from a restored iterate
  const int restore_iterations = 3;
  const auto data_dir = Test_Config::test_data_dir / "end-to-end_tests"
                        / "SingletScalarAllowed_test_nmax6";
  const auto dir_A = data_dir / "dual_feasible_jump" / "input" / "json";
  const auto dir_B = data_dir / "primal_feasible_jump" / "input" / "json";
  const std::string changed_A
    = "mat_huHER-bgnrovCSCBgKZ40XuFTcLhpW6hNmHx0gQEbwM.json";
  const std::string changed_B
    = "mat_RfFLGRJJdlsbUjbyx3vPfOJ2OVffJDUPiaJFv6v7euI.json";

  Test_Case_Runner runner("resident");
  const auto &output_dir = runner.output_dir;
  fs::create_directories(output_dir);
  Progress progress("resident");

  // PMP file lists: the changed block first, then the common files
  std::vector<std::string> common;
  for(const auto &entry : fs::directory_iterator(dir_A))
    if(entry.path().filename() != changed_A)
      common.push_back(fs::absolute(entry.path()).string());
  std::sort(common.begin(), common.end());
  REQUIRE(common.size() == 12);
  const auto write_nsv = [&](const fs::path &path, const fs::path &first) {
    std::string text = fs::absolute(first).string() + '\0';
    for(const auto &file : common)
      text += file + '\0';
    write_file(path, text);
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
      progress.step("pmp2sdp " + sdp.filename().string());
      runner.create_nested("pmp2sdp_" + sdp.filename().string())
        .mpi_run({Test_Config::program("pmp2sdp"), args}, num_procs);
    }

  // Block 0 (the changed block) is expensive: 2 ranks own it
  size_t num_blocks = 0;
  {
    const auto control_json = read_file(sdp_A / "control.json");
    const std::string key = "\"num_blocks\":";
    const auto pos = control_json.find(key);
    REQUIRE(pos != std::string::npos);
    num_blocks = std::stoul(control_json.substr(pos + key.size()));
  }
  REQUIRE(num_blocks == 11);
  const auto block_timings = block_timings_for(num_blocks);

  const auto sdpb_args_with = [&](const int max_iterations) {
    return boost::program_options::split_unix(
      "--checkpointInterval 3600 --maxRuntime 1341 "
      "--dualityGapThreshold 1.0e-30 --primalErrorThreshold 1.0e-200 "
      "--dualErrorThreshold 1.0e-200 --initialMatrixScalePrimal 1.0e20 "
      "--initialMatrixScaleDual 1.0e20 --feasibleCenteringParameter 0.1 "
      "--infeasibleCenteringParameter 0.3 --stepLengthReduction 0.7 "
      "--maxComplementarity 1.0e100 --maxIterations "
      + std::to_string(max_iterations)
      + " --verbosity 1 "
        "--writeSolution y,z --detectPrimalFeasibleJump "
        "--detectDualFeasibleJump "
        "--maxSharedMemory=100.1K --precision "
      + std::to_string(precision));
  };
  const auto sdpb_args = sdpb_args_with(1000);

  // File route: A, B, A through one checkpoint directory
  const auto ck_file = output_dir / "ck_file";
  fs::create_directories(ck_file);
  write_file(ck_file / "block_timings", block_timings);
  const auto ck_file_after_1 = output_dir / "ck_file_after_1";
  const std::vector<fs::path> sdps{sdp_A, sdp_B, sdp_A};
  for(size_t k = 0; k < sdps.size(); ++k)
    {
      const auto name = "file_" + std::to_string(k + 1);
      const Named_Args_Map args{
        {"--sdpDir", sdps[k].string()},
        {"--checkpointDir", ck_file.string()},
        {"--outDir", (output_dir / ("out_" + name)).string()}};
      progress.step("file route: sdpb, solve " + std::to_string(k + 1)
                    + " of 3 (" + (k == 1 ? "B" : "A") + ")");
      runner.create_nested("sdpb_" + name)
        .mpi_run({Test_Config::program("sdpb"), sdpb_args, args}, num_procs);
      snapshot_checkpoint(ck_file, output_dir / ("ck_" + name));
      if(k == 0)
        fs::copy(ck_file, ck_file_after_1, fs::copy_options::recursive);
    }
  // B from the checkpoint of the first solve, stopped after
  // restore_iterations
  {
    const Named_Args_Map args{
      {"--sdpDir", sdp_B.string()},
      {"--checkpointDir", ck_file_after_1.string()},
      {"--outDir", (output_dir / "out_file_restored").string()}};
    progress.step("file route: sdpb, B from the first solve's checkpoint, "
                  + std::to_string(restore_iterations) + " iterations");
    runner.create_nested("sdpb_file_restored")
      .mpi_run({Test_Config::program("sdpb"),
                sdpb_args_with(restore_iterations), args},
               num_procs);
    snapshot_checkpoint(ck_file_after_1, output_dir / "ck_file_restored");
  }

  // Resident: the same solves, commands written in advance
  const auto ck_res = output_dir / "ck_res";
  fs::create_directories(ck_res);
  write_file(ck_res / "block_timings", block_timings);
  const auto control = output_dir / "control";
  fs::create_directories(control);
  const auto solve = [&](int k, const std::string &extra_args = "") {
    return "{\"cmd\": \"solve\", \"checkpoint\": true, \"args\": [\"--outDir="
           + (output_dir / ("out_res_" + std::to_string(k))).string()
           + "\", \"--checkpointDir="
           + (output_dir / ("ck_res_" + std::to_string(k))).string() + "\""
           + extra_args + "]}";
  };
  const auto ck_saved = output_dir / "ck_saved_1";
  write_file(control / "cmd.1.json",
             "{\"commands\": [" + solve(1)
               + ", {\"cmd\": \"save_iterate\", \"name\": \"after_1\"}, "
                 "{\"cmd\": \"save_iterate\", \"dir\": \""
               + ck_saved.string() + "\"}]}");
  write_file(control / "cmd.2.json",
             "{\"commands\": [{\"cmd\": \"update_pmp\", \"pmp\": \""
               + nsv_B.string() + "\"}, " + solve(2) + "]}");
  write_file(control / "cmd.3.json",
             "{\"commands\": [{\"cmd\": \"update_sdp\", \"sdp\": \""
               + sdp_A.string() + "\"}, " + solve(3) + "]}");
  // "inherit": false takes nothing from the start-up command line: without
  // its own --precision the solve is refused (with inheritance it would
  // run at the start-up precision)
  write_file(control / "cmd.4.json",
             "{\"cmd\": \"solve\", \"inherit\": false, \"args\": "
             "[\"--sdpDir="
               + sdp_A.string() + "\", \"--outDir="
               + (output_dir / "out_res_inherit").string() + "\"]}");
  // B again, from the iterate saved after solve 1, three ways, each for
  // restore_iterations iterations
  const std::string restore_args
    = ", \"--maxIterations=" + std::to_string(restore_iterations) + "\"";
  const auto solve_from = [&](int k, const std::string &start) {
    auto text = solve(k, restore_args);
    text.insert(text.size() - 1, ", \"start\": " + start);
    return text;
  };
  write_file(control / "cmd.5.json",
             "{\"commands\": [{\"cmd\": \"update_sdp\", \"sdp\": \""
               + sdp_B.string()
               + "\"}, {\"cmd\": \"restore_iterate\", \"name\": "
                 "\"after_1\"}, "
               + solve(4, restore_args) + "]}");
  write_file(control / "cmd.6.json",
             "{\"commands\": [{\"cmd\": \"restore_iterate\", \"dir\": \""
               + ck_saved.string() + "\"}, " + solve(5, restore_args) + "]}");
  write_file(control / "cmd.7.json",
             "{\"commands\": ["
               + solve_from(6, "{\"iterate\": \"after_1\"}") + "]}");
  write_file(control / "cmd.8.json", "{\"cmd\": \"quit\"}");
  const Named_Args_Map resident_args{{"--sdpDir", sdp_A.string()},
                                     {"--checkpointDir", ck_res.string()},
                                     {"--controlDir", control.string()}};
  progress.step("resident: sdpb --controlDir, 8 command files");
  progress.watch(control);
  runner.create_nested("sdpb_resident")
    .mpi_run({Test_Config::program("sdpb"), sdpb_args, resident_args}, num_procs);
  progress.step("comparing the outputs");
  REQUIRE(fs::exists(control / "exit.json"));
  REQUIRE(read_file(control / "ready.json").find("\"solve_inherit\"")
          != std::string::npos);
  for(int seq = 1; seq <= 8; ++seq)
    {
      CAPTURE(seq);
      const auto result
        = read_file(control / ("res." + std::to_string(seq) + ".json"));
      CAPTURE(result);
      if(seq == 4)
        {
          REQUIRE(result.find("\"error\"") != std::string::npos);
          REQUIRE(result.find("--precision is fixed") != std::string::npos);
        }
      else
        REQUIRE(result.find("\"error\"") == std::string::npos);
    }

  // The solves switch between excluded and allowed
  REQUIRE(terminate_reason(output_dir / "out_file_1")
          != terminate_reason(output_dir / "out_file_2"));
  REQUIRE(terminate_reason(output_dir / "out_file_1")
          == terminate_reason(output_dir / "out_file_3"));

  for(int k = 1; k <= 3; ++k)
    {
      CAPTURE(k);
      const auto out_file = output_dir / ("out_file_" + std::to_string(k));
      const auto out_res = output_dir / ("out_res_" + std::to_string(k));
      REQUIRE(out_txt_without_runtime(out_file)
              == out_txt_without_runtime(out_res));
      for(const auto *name : {"y.txt", "z.txt"})
        {
          CAPTURE(name);
          REQUIRE(read_file(out_file / name) == read_file(out_res / name));
        }
      // Generation k in both: one checkpoint per solve
      const auto ck_f = output_dir / ("ck_file_" + std::to_string(k));
      const auto ck_r = output_dir / ("ck_res_" + std::to_string(k));
      const auto files = checkpoint_files(ck_f, k);
      REQUIRE(files.size() == num_procs);
      REQUIRE(checkpoint_files(ck_r, k) == files);
      for(const auto &name : files)
        {
          CAPTURE(name);
          REQUIRE(read_file(ck_f / name) == read_file(ck_r / name));
        }
    }

  // Solves 4, 5, 6 (B from the iterate saved after solve 1) equal the file
  // route's solve of B from the checkpoint of solve 1, all stopped after
  // restore_iterations. The iterate after solve 3 differs from the saved one
  // on every rank, so a restore that kept part of it would differ.
  {
    const auto saved = checkpoint_files_by_rank(output_dir / "ck_file_1", 1);
    const auto third = checkpoint_files_by_rank(output_dir / "ck_res_3", 3);
    REQUIRE(saved.size() == num_procs);
    REQUIRE(third.size() == num_procs);
    for(const auto &[rank_name, path] : saved)
      {
        CAPTURE(rank_name);
        REQUIRE(read_file(path) != read_file(third.at(rank_name)));
      }
  }
  const auto out_restored = output_dir / "out_file_restored";
  REQUIRE(terminate_reason(out_restored).find("maxIterations exceeded")
          != std::string::npos);
  const auto ck_f2
    = checkpoint_files_by_rank(output_dir / "ck_file_restored", 2);
  REQUIRE(ck_f2.size() == num_procs);
  for(int k = 4; k <= 6; ++k)
    {
      CAPTURE(k);
      const auto &out_file = out_restored;
      const auto out_res = output_dir / ("out_res_" + std::to_string(k));
      REQUIRE(out_txt_without_runtime(out_file)
              == out_txt_without_runtime(out_res));
      for(const auto *name : {"y.txt", "z.txt"})
        {
          CAPTURE(name);
          REQUIRE(read_file(out_file / name) == read_file(out_res / name));
        }
      const auto ck_r
        = checkpoint_files_by_rank(output_dir / ("ck_res_" + std::to_string(k)));
      REQUIRE(ck_r.size() == ck_f2.size());
      for(const auto &[rank_name, path] : ck_f2)
        {
          CAPTURE(rank_name);
          REQUIRE(ck_r.count(rank_name) == 1);
          REQUIRE(read_file(path) == read_file(ck_r.at(rank_name)));
        }
    }
}

TEST_CASE("resident_family")
{
  INFO("sdpb --controlDir: a parametric block (load_family, "
       "set_coefficients) with an exact contraction equals the file route");
  // Block 0 is lambda^T G lambda >= 0, a 1x1 constant constraint with a 2x2 G:
  // G_00 = unit (the normalization), G_01 = 0, G_11 = -unit. With unit
  // lambdas, (1,0) and (0,1), the block is G_00 or G_11 exactly, so the
  // resident's contraction must equal pmp2sdp's conversion of that component.
  const int num_procs = 4;
  const int precision = 384;
  const auto dir_A = Test_Config::test_data_dir / "end-to-end_tests"
                     / "SingletScalarAllowed_test_nmax6" / "dual_feasible_jump"
                     / "input" / "json";
  const std::string norm_file
    = "mat__jiutIpgPLv34k9-62shY0ddvi3BmZGH3YmLnKhtYAE.json";

  Test_Case_Runner runner("resident_family");
  const auto &output_dir = runner.output_dir;
  fs::create_directories(output_dir);
  Progress progress("resident_family");

  // The normalization strings, read as text
  std::vector<std::string> unit;
  {
    rapidjson::Document doc;
    doc.Parse(read_file(dir_A / norm_file).c_str());
    REQUIRE(doc.HasMember("normalization"));
    for(const auto &v : doc["normalization"].GetArray())
      unit.emplace_back(v.GetString());
  }
  const auto const_pmp = [](const std::vector<std::string> &vec) {
    std::string text = "{\"PositiveMatrixWithPrefactorArray\": [{\"prefactor\": "
                       "{\"base\": \"1\", \"constant\": \"1\", \"poles\": []}, "
                       "\"polynomials\": [[[";
    for(size_t i = 0; i < vec.size(); ++i)
      text += (i == 0 ? "[\"" : ", [\"") + vec[i] + "\"]";
    return text + "]]]}]}";
  };
  std::vector<std::string> minus_unit, zero(unit.size(), "0");
  for(const auto &x : unit)
    minus_unit.push_back(x.rfind('-', 0) == 0 ? x.substr(1) : "-" + x);
  const auto comp_00 = fs::absolute(output_dir / "comp_00.json");
  const auto comp_01 = fs::absolute(output_dir / "comp_01.json");
  const auto comp_11 = fs::absolute(output_dir / "comp_11.json");
  write_file(comp_00, const_pmp(unit));
  write_file(comp_01, const_pmp(zero));
  write_file(comp_11, const_pmp(minus_unit));

  // The base problem (all its files), with the family block first
  std::vector<std::string> common;
  for(const auto &entry : fs::directory_iterator(dir_A))
    common.push_back(fs::absolute(entry.path()).string());
  std::sort(common.begin(), common.end());
  const auto write_nsv = [&](const fs::path &path, const fs::path &first) {
    std::string text = first.string() + '\0';
    for(const auto &file : common)
      text += file + '\0';
    write_file(path, text);
  };
  const auto nsv_X = output_dir / "pmp_X.nsv";
  const auto nsv_Y = output_dir / "pmp_Y.nsv";
  write_nsv(nsv_X, comp_00);
  write_nsv(nsv_Y, comp_11);
  const auto sdp_X = output_dir / "sdp_X";
  const auto sdp_Y = output_dir / "sdp_Y";
  for(const auto &[nsv, sdp] :
      {std::pair{nsv_X, sdp_X}, std::pair{nsv_Y, sdp_Y}})
    {
      const Named_Args_Map args{{"--input", nsv.string()},
                                {"--output", sdp.string()},
                                {"--precision", std::to_string(precision)}};
      progress.step("pmp2sdp " + sdp.filename().string());
      runner.create_nested("pmp2sdp_" + sdp.filename().string())
        .mpi_run({Test_Config::program("pmp2sdp"), args}, num_procs);
    }
  size_t num_blocks = 0;
  {
    const auto control_json = read_file(sdp_X / "control.json");
    const std::string key = "\"num_blocks\":";
    const auto pos = control_json.find(key);
    REQUIRE(pos != std::string::npos);
    num_blocks = std::stoul(control_json.substr(pos + key.size()));
  }
  const auto block_timings = block_timings_for(num_blocks);

  const auto sdpb_args = boost::program_options::split_unix(
    "--checkpointInterval 3600 --maxRuntime 1341 "
    "--dualityGapThreshold 1.0e-30 --primalErrorThreshold 1.0e-200 "
    "--dualErrorThreshold 1.0e-200 --initialMatrixScalePrimal 1.0e20 "
    "--initialMatrixScaleDual 1.0e20 --feasibleCenteringParameter 0.1 "
    "--infeasibleCenteringParameter 0.3 --stepLengthReduction 0.7 "
    "--maxComplementarity 1.0e100 --maxIterations 1000 --verbosity 1 "
    "--writeSolution y,z --detectPrimalFeasibleJump --detectDualFeasibleJump "
    "--precision "
    + std::to_string(precision));

  // File route: X, Y, X
  const auto ck_file = output_dir / "ck_file";
  fs::create_directories(ck_file);
  write_file(ck_file / "block_timings", block_timings);
  const std::vector<fs::path> sdps{sdp_X, sdp_Y, sdp_X};
  for(size_t k = 0; k < sdps.size(); ++k)
    {
      const auto name = "file_" + std::to_string(k + 1);
      const Named_Args_Map args{
        {"--sdpDir", sdps[k].string()},
        {"--checkpointDir", ck_file.string()},
        {"--outDir", (output_dir / ("out_" + name)).string()}};
      progress.step("file route: sdpb, solve " + std::to_string(k + 1)
                    + " of 3 (" + (k == 1 ? "Y" : "X") + ")");
      runner.create_nested("sdpb_" + name)
        .mpi_run({Test_Config::program("sdpb"), sdpb_args, args}, num_procs);
      snapshot_checkpoint(ck_file, output_dir / ("ck_" + name));
    }

  // Resident: the family, lambda = (1,0), (0,1), (1,0)
  const auto ck_res = output_dir / "ck_res";
  fs::create_directories(ck_res);
  write_file(ck_res / "block_timings", block_timings);
  const auto control = output_dir / "control";
  fs::create_directories(control);
  const auto solve = [&](int k) {
    return "{\"cmd\": \"solve\", \"checkpoint\": true, \"args\": [\"--outDir="
           + (output_dir / ("out_res_" + std::to_string(k))).string()
           + "\", \"--checkpointDir="
           + (output_dir / ("ck_res_" + std::to_string(k))).string() + "\"]}";
  };
  const auto set = [](const std::string &l0, const std::string &l1) {
    return "{\"cmd\": \"set_coefficients\", \"family\": \"ext\", "
           "\"lambda\": [\""
           + l0 + "\", \"" + l1 + "\"]}";
  };
  write_file(control / "cmd.1.json",
             "{\"commands\": [{\"cmd\": \"load_family\", \"name\": \"ext\", "
             "\"block\": 0, \"form_dim\": 2, \"normalization_file\": \""
               + fs::absolute(dir_A / norm_file).string()
               + "\", \"components\": [\"" + comp_00.string() + "\", \""
               + comp_01.string() + "\", \"" + comp_11.string() + "\"]}, "
               + set("1", "0") + ", " + solve(1) + "]}");
  write_file(control / "cmd.2.json", "{\"commands\": [" + set("0", "1") + ", "
                                       + solve(2) + "]}");
  write_file(control / "cmd.3.json", "{\"commands\": [" + set("1", "0") + ", "
                                       + solve(3) + "]}");
  write_file(control / "cmd.4.json", "{\"cmd\": \"quit\"}");
  const Named_Args_Map resident_args{{"--sdpDir", sdp_X.string()},
                                     {"--checkpointDir", ck_res.string()},
                                     {"--controlDir", control.string()}};
  progress.step("resident: sdpb --controlDir, 4 command files");
  progress.watch(control);
  runner.create_nested("sdpb_resident")
    .mpi_run({Test_Config::program("sdpb"), sdpb_args, resident_args}, num_procs);
  progress.step("comparing the outputs");
  for(int seq = 1; seq <= 4; ++seq)
    {
      CAPTURE(seq);
      const auto result
        = read_file(control / ("res." + std::to_string(seq) + ".json"));
      CAPTURE(result);
      REQUIRE(result.find("\"error\"") == std::string::npos);
    }
  REQUIRE(terminate_reason(output_dir / "out_file_1")
          != terminate_reason(output_dir / "out_file_2"));

  for(int k = 1; k <= 3; ++k)
    {
      CAPTURE(k);
      const auto out_file = output_dir / ("out_file_" + std::to_string(k));
      const auto out_res = output_dir / ("out_res_" + std::to_string(k));
      REQUIRE(out_txt_without_runtime(out_file)
              == out_txt_without_runtime(out_res));
      for(const auto *name : {"y.txt", "z.txt"})
        {
          CAPTURE(name);
          REQUIRE(read_file(out_file / name) == read_file(out_res / name));
        }
      const auto ck_f = output_dir / ("ck_file_" + std::to_string(k));
      const auto ck_r = output_dir / ("ck_res_" + std::to_string(k));
      const auto files = checkpoint_files(ck_f, k);
      REQUIRE(files.size() == num_procs);
      REQUIRE(checkpoint_files(ck_r, k) == files);
      for(const auto &name : files)
        {
          CAPTURE(name);
          REQUIRE(read_file(ck_f / name) == read_file(ck_r / name));
        }
    }
}
