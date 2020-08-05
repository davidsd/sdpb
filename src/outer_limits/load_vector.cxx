#include <El.hpp>
#include <boost/filesystem.hpp>
#include <boost/filesystem/fstream.hpp>

std::vector<El::BigFloat>
load_vector(const boost::filesystem::path &vector_path)
{
  boost::filesystem::ifstream vector_file(vector_path);
  std::string vector_string;
  vector_file >> vector_string;
  std::vector<El::BigFloat> result;
  while(vector_file)
    {
      result.emplace_back(vector_string);
      vector_file >> vector_string;
    }
  return result;
}
