#pragma once

#include <catch2/catch_amalgamated.hpp>
#include <boost/filesystem.hpp>
#include <string>

namespace Test_Config
{
  inline std::string mpirun = "mpirun";
  const boost::filesystem::path test_data_dir = "test/data";
  const boost::filesystem::path test_output_dir = "test/out";
  const boost::filesystem::path test_log_dir = "test/out/log";
}
