#include "Function_Block.hxx"

#include "../../sdp_read.hxx"

#include <boost/filesystem.hpp>

void read_json(const boost::filesystem::path &input_path,
               std::vector<El::BigFloat> &objectives,
               std::vector<El::BigFloat> &normalization,
               std::vector<std::vector<std::map<El::BigFloat, El::BigFloat>>> &functions);

void read_function_blocks(
  const boost::filesystem::path &input_file,
  std::vector<El::BigFloat> &objectives,
  std::vector<El::BigFloat> &normalization,
  std::vector<std::vector<std::map<El::BigFloat, El::BigFloat>>> &functions)
{
  if(input_file.extension() == ".nsv")
    {
      for(auto &filename : read_file_list(input_file))
        {
          if(!filename.empty())
            {
              read_function_blocks(filename, objectives, normalization,
                                   functions);
            }
        }
    }
  else if(input_file.extension() == ".json")
    {
      read_json(input_file, objectives, normalization, functions);
    }
}
