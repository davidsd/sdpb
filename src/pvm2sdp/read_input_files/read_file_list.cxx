#include <boost/filesystem.hpp>
#include <boost/filesystem/fstream.hpp>
#include <boost/algorithm/string.hpp>

std::vector<boost::filesystem::path>
read_file_list(const boost::filesystem::path &filename)
{
  auto file_size(boost::filesystem::file_size(filename));
  std::vector<char> file_contents(file_size);
  file_contents.reserve(file_size);
  boost::filesystem::ifstream infile(filename);
  infile.read(file_contents.data(),file_size);
  std::vector<boost::filesystem::path> result;
  boost::split(result,file_contents,boost::is_any_of("\0"));
  return result;
}

