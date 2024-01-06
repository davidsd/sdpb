#include "PMWP_File_Parse_Result.hxx"

#include "read_json/read_json.hxx"
#include "read_mathematica/read_mathematica.hxx"

namespace fs = std::filesystem;

PMWP_File_Parse_Result::PMWP_File_Parse_Result(
  const fs::path &input_path,
  const std::function<bool(size_t matrix_index)> &should_parse_matrix)
    : path(input_path), num_matrices(0)
{
  if(input_path.extension() == ".json")
    {
      read_json(input_path, should_parse_matrix, objective, normalization,
                num_matrices, parsed_matrices);
    }
  else if(input_path.extension() == ".m")
    {
      read_mathematica(input_path, should_parse_matrix, objective,
                       normalization, num_matrices, parsed_matrices);
    }
  else
    {
      El::RuntimeError("Cannot parse input file: ", input_path,
                       ". Expected .json or .m extension.");
    }
  for(auto &[index, matrix] : parsed_matrices)
    {
      if(index >= num_matrices)
        El::LogicError("index=", index, "should be less than", num_matrices,
                       ", input_path=", input_path);
    }
}