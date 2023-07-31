#include <boost/filesystem.hpp>
#include <boost/filesystem/fstream.hpp>
#include <boost/algorithm/string.hpp>

#include <stdexcept>
#include <vector>

// Read .nsv containing null-separated list of files
std::vector<boost::filesystem::path>
read_nsv_file_list(const boost::filesystem::path &input_file)
{
  auto file_size(boost::filesystem::file_size(input_file));
  std::vector<char> file_contents(file_size);
  file_contents.reserve(file_size);
  boost::filesystem::ifstream infile(input_file);
  infile.read(file_contents.data(), file_size);
  if(!infile.good())
    {
      throw std::runtime_error("Unable to read: " + input_file.string());
    }
  std::vector<boost::filesystem::path> result;
  using namespace std::literals;
  boost::split(result, file_contents, boost::is_any_of("\0"s));
  result.erase(std::remove_if(result.begin(), result.end(),
                              [](const boost::filesystem::path &path) {
                                return path.empty();
                              }),
               result.end());
  for(auto &path : result)
    {
      if(path.is_relative())
        path = input_file.parent_path() / path;
    }
  return result;
}
