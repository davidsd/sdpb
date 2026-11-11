#pragma once

#include "process.hxx"

#include <catch2/catch_amalgamated.hpp>
#include <filesystem>
#include <string>

namespace Test_Config
{
  inline Command mpirun{"mpirun", {}, {}};
  // Directory with SDPB executables, set in main()
  inline std::filesystem::path build_dir;
  // Path to SDPB executable, e.g. program("sdpb") = "/path/to/build/sdpb"
  inline std::string program(const std::string &name)
  {
    return (build_dir / name).string();
  }
  const std::filesystem::path test_data_dir = "test/data";
  // Output is in the build directory, so that different builds do not overwrite it
  inline std::filesystem::path test_output_dir()
  {
    return build_dir / "test" / "out";
  }
  inline std::filesystem::path test_log_dir()
  {
    return test_output_dir() / "log";
  }
}
