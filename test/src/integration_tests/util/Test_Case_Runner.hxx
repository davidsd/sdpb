#pragma once

#include "process.hxx"

#include <catch2/catch_amalgamated.hpp>

#include <filesystem>
#include <string>
#include <boost/noncopyable.hpp>

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
    using Named_Args_Map = Command::Named_Args_Map;

    const std::string name;
    const std::filesystem::path data_dir;
    const std::filesystem::path output_dir;
    const std::filesystem::path stdout_path;
    const std::filesystem::path stderr_path;

    explicit Test_Case_Runner(const std::string &name);

    // create Test_Case_Runner "name/suffix"
    [[nodiscard]] Test_Case_Runner
    create_nested(const std::string &suffix,
                  const std::string &separator = "/") const;

    void run(const Command &command, int required_exit_code = 0,
             const std::string &required_error_msg = "") const;
    void mpi_run(const Command &command, int numProcs = 2,
                 int required_exit_code = 0,
                 const std::string &required_error_msg = "") const;

    [[nodiscard]] std::filesystem::path
    unzip_to_temp_dir(const std::filesystem::path &zip_path) const;
  };
}
