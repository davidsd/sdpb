#include "JSON_Parser.hxx"

#include "../../../ostream_vector.hxx"

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

  std::cout << parser.objective_state.value << "\n";
  std::cout << parser.normalization_state.value << "\n";
  std::cout << parser.positive_matrices_with_prefactor_state.element_state.damped_rational_state.value.constant << "\n";
  std::cout << parser.positive_matrices_with_prefactor_state.element_state.damped_rational_state.value.base << "\n";
  std::cout << parser.positive_matrices_with_prefactor_state.element_state.damped_rational_state.value.poles << "\n";
  exit(0);
}

