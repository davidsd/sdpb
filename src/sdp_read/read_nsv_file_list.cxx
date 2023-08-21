#include <filesystem>
#include <boost/algorithm/string.hpp>

#include <stdexcept>
#include <vector>
#include <fstream>

namespace fs = std::filesystem;

// Read .nsv containing null-separated list of files
std::vector<fs::path> read_nsv_file_list(const fs::path &input_file)
{
  auto file_size(fs::file_size(input_file));
  std::vector<char> file_contents(file_size);
  file_contents.reserve(file_size);
  std::ifstream infile(input_file);
  infile.read(file_contents.data(), file_size);
  if(!infile.good())
    {
      throw std::runtime_error("Unable to read: " + input_file.string());
    }
  std::vector<fs::path> result;
  using namespace std::literals;
  boost::split(result, file_contents, boost::is_any_of("\0"s));
  result.erase(std::remove_if(result.begin(), result.end(),
                   [](const fs::path &path) { return path.empty();
                              }),
               result.end());
  for(auto &path : result)
    {
      if(path.is_relative())
        path = input_file.parent_path() / path;
    }
  return result;
}
