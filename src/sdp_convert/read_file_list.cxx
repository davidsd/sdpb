#include <boost/filesystem.hpp>
#include <boost/filesystem/fstream.hpp>
#include <boost/algorithm/string.hpp>

#include <iostream>
#include <stdexcept>
#include <vector>

std::vector<boost::filesystem::path>
read_file_list(const boost::filesystem::path &filename)
{
  auto file_size(boost::filesystem::file_size(filename));
  std::vector<char> file_contents(file_size);
  file_contents.reserve(file_size);
  boost::filesystem::ifstream infile(filename);
  infile.read(file_contents.data(), file_size);
  if(!infile.good())
    {
      throw std::runtime_error("Unable to read: " + filename.string());
    }
  std::vector<boost::filesystem::path> result;
  using namespace std::literals;
  boost::split(result, file_contents, boost::is_any_of("\0"s));
  return result;
}
