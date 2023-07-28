#pragma once

#include <catch2/catch_amalgamated.hpp>

#include <boost/filesystem.hpp>
#include <boost/process.hpp>
#include <sstream>
#include <string>
#include <variant>

namespace Test_Config
{
  inline std::string mpirun = "mpirun";
  const boost::filesystem::path test_data_dir = "test/data";
  const boost::filesystem::path test_output_dir = "test/out";
  const boost::filesystem::path test_log_dir = "test/out/log";
}

namespace Test_Util
{
  typedef std::map<std::string, std::string> Named_Args_Map;

  // Sets up and cleans working directories.
  // Runs commands (see run() and mpi_run() methods) and redirects output to
  // "$test_log_dir/$name.stdout.log"
  //
  // NB: if TEST_CASE has several SECTION's, create Test_Case_Runner only
  // inside the section! Otherwise, some output files will be removed.
  struct Test_Case_Runner : boost::noncopyable
  {
    const std::string name;
    const boost::filesystem::path data_dir;
    const boost::filesystem::path output_dir;
    const boost::filesystem::path stdout_path;
    const boost::filesystem::path stderr_path;

    explicit Test_Case_Runner(const std::string &name);

    // create Test_Case_Runner "name/suffix"
    Test_Case_Runner create_nested(const std::string &suffix,
                                   const std::string &separator = "/") const;

    [[nodiscard]] int run(const std::string &command) const;
    [[nodiscard]] int
    mpi_run(const std::string &command, int numProcs = 2) const;
    [[nodiscard]] int run(const std::vector<std::string> &args,
                          const Named_Args_Map &named_args = {}) const;
    [[nodiscard]] int
    mpi_run(const std::vector<std::string> &args,
            const Named_Args_Map &named_args = {}, int numProcs = 2) const;

    [[nodiscard]] int diff(const boost::filesystem::path &x,
                           const boost::filesystem::path &y) const;

    // diff two zip archives ignoring control.json file
    // control.json may differ by "command" field - we ignore it.
    [[nodiscard]] int diff_sdp_zip(const boost::filesystem::path &x,
                                 const boost::filesystem::path &y) const;

    // Check if stderr contains given substring.
    // NB: multiline strings not supported!
    [[nodiscard]] bool
    stderr_contains_substring(const std::string &substring) const;
  };
}
