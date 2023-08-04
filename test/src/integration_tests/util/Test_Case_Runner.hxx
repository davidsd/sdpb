#pragma once

#include <catch2/catch_amalgamated.hpp>

#include <boost/filesystem.hpp>
#include <boost/process.hpp>
#include <string>

namespace Test_Util
{
  // Sets up and cleans working directories.
  // Runs commands (see run() and mpi_run() methods) and redirects output to
  // "$test_log_dir/$name.stdout.log"
  //
  // NB: if TEST_CASE has several SECTION's, create Test_Case_Runner only
  // inside the section! Otherwise, some output files will be removed.
  struct Test_Case_Runner : boost::noncopyable
  {
    typedef std::map<std::string, std::string> Named_Args_Map;

    const std::string name;
    const boost::filesystem::path data_dir;
    const boost::filesystem::path output_dir;
    const boost::filesystem::path stdout_path;
    const boost::filesystem::path stderr_path;

    explicit Test_Case_Runner(const std::string &name);

    // create Test_Case_Runner "name/suffix"
    [[nodiscard]] Test_Case_Runner
    create_nested(const std::string &suffix,
                  const std::string &separator = "/") const;

    [[nodiscard]] int run(const std::string &command) const;
    [[nodiscard]] int
    mpi_run(const std::string &command, int numProcs = 2) const;
    [[nodiscard]] int run(const std::vector<std::string> &args,
                          const Named_Args_Map &named_args = {}) const;
    [[nodiscard]] int
    mpi_run(const std::vector<std::string> &args,
            const Named_Args_Map &named_args = {}, int numProcs = 2) const;

    [[nodiscard]] boost::filesystem::path
    unzip_to_temp_dir(const boost::filesystem::path &zip_path) const;

    // Check if stderr contains given substring.
    // NB: multiline strings not supported!
    [[nodiscard]] bool
    stderr_contains_substring(const std::string &substring) const;
  };
}
