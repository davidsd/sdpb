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
  reader.Parse(wrapper, parser);

  if(!parser.objective_state.value.empty())
    {
      std::swap(objectives, parser.objective_state.value);
    }
  if(!parser.normalization_state.value.empty())
    {
      std::swap(normalization, parser.normalization_state.value);
    }
  size_t offset(matrices.size());
  auto &temp_matrices(
    parser.positive_matrices_with_prefactor_state.value);
  matrices.resize(matrices.size() + temp_matrices.size());
  
  size_t rank(El::mpi::Rank(El::mpi::COMM_WORLD)),
    num_procs(El::mpi::Size(El::mpi::COMM_WORLD));
  for(size_t index = 0; index < temp_matrices.size(); ++index)
    {
      if((offset + index) % num_procs == rank)
        {
          std::swap(matrices[offset + index], temp_matrices[index]);
        }
    }
}
