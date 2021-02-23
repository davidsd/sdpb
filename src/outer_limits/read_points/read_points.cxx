#include "../../sdp_read.hxx"

#include <boost/filesystem.hpp>

void read_points_json(const boost::filesystem::path &input_path,
                      std::vector<std::vector<El::BigFloat>> &points);

void read_points(const boost::filesystem::path &input_path,
                 std::vector<std::vector<El::BigFloat>> &points)
{
  if(input_path.extension() == ".nsv")
    {
      for(auto &filename : read_file_list(input_path))
        {
          read_points(filename, points);
        }
    }
  else if(input_path.extension() == ".json")
    {
      read_points_json(input_path, points);
    }
  else
    {
      throw std::runtime_error("Unknown extension for file: "
                               + input_path.string());
    }
}
