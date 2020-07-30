#include <El.hpp>
#include <boost/filesystem.hpp>
#include <boost/filesystem/fstream.hpp>

std::vector<El::BigFloat>
load_objective(const boost::filesystem::path &objective_path)
{
  boost::filesystem::ifstream objective_file(objective_path);
  std::string objective_string;
  objective_file >> objective_string;
  std::vector<El::BigFloat> result;
  while(objective_file)
    {
      result.emplace_back(objective_string);
      objective_file >> objective_string;
    }
  return result;
}
