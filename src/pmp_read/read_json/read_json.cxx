#include "read_json.hxx"
#include "JSON_Parser.hxx"
#include "pmp/Polynomial_Vector_Matrix.hxx"
#include "sdpb_util/to_matrix.hxx"

#include <rapidjson/istreamwrapper.h>

namespace fs = std::filesystem;

void read_json(
  const std::filesystem::path &input_path,
  const std::function<bool(size_t matrix_index)> &should_parse_matrix,
  std::vector<El::BigFloat> &objectives,
  std::vector<El::BigFloat> &normalization, size_t &num_matrices,
  std::map<size_t, Polynomial_Vector_Matrix> &parsed_matrices)
{
  std::ifstream input_file(input_path);
  rapidjson::IStreamWrapper wrapper(input_file);
  JSON_Parser parser;
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

  // TODO skip matrices that we don't need
  auto &all_matrices(parser.positive_matrices_with_prefactor_state.value);
  num_matrices = all_matrices.size();

  for(size_t index = 0; index < all_matrices.size(); ++index)
    {
      if(should_parse_matrix(index))
        {
          parsed_matrices.emplace(index, std::move(*all_matrices.at(index)));
        }
    }
}
