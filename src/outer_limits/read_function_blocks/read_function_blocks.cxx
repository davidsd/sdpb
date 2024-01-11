#include "pmp_read/pmp_read.hxx"
#include "outer_limits/Function.hxx"

#include <filesystem>

namespace fs = std::filesystem;

void read_json(
  const fs::path &input_path, std::vector<El::BigFloat> &objectives,
  std::vector<El::BigFloat> &normalization,
  std::vector<std::vector<std::vector<std::vector<Function>>>> &functions);

void read_function_blocks(
  const fs::path &input_file, std::vector<El::BigFloat> &objectives,
  std::vector<El::BigFloat> &normalization,
  std::vector<std::vector<std::vector<std::vector<Function>>>> &functions)
{
  if(input_file.extension() == ".nsv")
    {
      for(auto &filename : read_nsv_file_list(input_file))
        {
          read_function_blocks(filename, objectives, normalization, functions);
        }
    }
  else if(input_file.extension() == ".json")
    {
      read_json(input_file, objectives, normalization, functions);
    }
}
