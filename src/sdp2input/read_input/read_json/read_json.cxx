#include "JSON_Parser.hxx"

#include <rapidjson/istreamwrapper.h>
#include <boost/filesystem/fstream.hpp>

void read_json(const boost::filesystem::path &input_path,
               std::vector<El::BigFloat> &objectives,
               std::vector<El::BigFloat> &normalization,
               std::vector<Positive_Matrix_With_Prefactor> &matrices)
{
  boost::filesystem::ifstream input_file(input_path);
  rapidjson::IStreamWrapper wrapper(input_file);
  JSON_Parser parser;
  rapidjson::Reader reader;
  reader.Parse(wrapper,parser);
  exit(0);
}

