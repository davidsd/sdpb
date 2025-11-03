#pragma once

#include "process.hxx"

#include <catch2/catch_amalgamated.hpp>
#include <filesystem>
#include <string>

namespace Test_Config
{
  inline Command mpirun{"mpirun", {}, {}};
  const std::filesystem::path test_data_dir = "test/data";
  const std::filesystem::path test_output_dir = "test/out";
  const std::filesystem::path test_log_dir = "test/out/log";
}
