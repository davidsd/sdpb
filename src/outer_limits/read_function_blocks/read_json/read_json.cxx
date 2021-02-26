#include "Function_Blocks_Parser.hxx"

#include <rapidjson/istreamwrapper.h>
#include <boost/filesystem.hpp>
#include <boost/filesystem/fstream.hpp>

void read_json(
  const boost::filesystem::path &input_path,
  std::vector<El::BigFloat> &objectives,
  std::vector<El::BigFloat> &normalization,
  std::vector<std::vector<std::vector<std::vector<Function>>>> &functions)
{
  boost::filesystem::ifstream input_file(input_path);
  rapidjson::IStreamWrapper wrapper(input_file);
  Function_Blocks_Parser parser;
  rapidjson::Reader reader;
  reader.Parse(wrapper, parser);

  if(!parser.objective_state.value.empty())
    {
      std::swap(objectives, parser.objective_state.value);
    }
  if(!parser.normalization_state.value.empty())
    {
      std::swap(normalization, parser.normalization_state.value);
    }

  size_t functions_offset(functions.size());
  auto &temp_functions(parser.functions_state.value);
  functions.resize(functions.size() + temp_functions.size());
  for(size_t index(0); index != temp_functions.size(); ++index)
    {
      std::swap(functions[functions_offset + index], temp_functions[index]);
    }
}
