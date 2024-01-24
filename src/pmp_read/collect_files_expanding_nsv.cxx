#include "pmp_read.hxx"

#include <filesystem>
#include <variant>
#include <vector>

namespace fs = std::filesystem;

// collect_files_expanding_nsv() replaces
// each .nsv in input_files with its content.
// .nsv inside .nsv are also expanded.

namespace
{
  void collect_files_expanding_nsv(const fs::path &input_file,
                                   std::vector<fs::path> &output)
  {
    if(input_file.empty())
      return;
    if(input_file.extension() == ".nsv")
      {
        for(const auto &inner_file : read_nsv_file_list(input_file))
          {
            collect_files_expanding_nsv(inner_file, output);
          }
      }
    else
      {
        output.push_back(input_file);
      }
  }
}

std::vector<fs::path>
collect_files_expanding_nsv(const std::vector<fs::path> &input_files)
{
  std::vector<fs::path> output;
  for(const auto &file : input_files)
    {
      collect_files_expanding_nsv(file, output);
    }
  return output;
}

std::vector<fs::path> collect_files_expanding_nsv(const fs::path &input_file)
{
  std::vector<fs::path> output;
  collect_files_expanding_nsv(input_file, output);
  return output;
}
