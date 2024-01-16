#include "PMP_File_Parse_Result.hxx"

namespace fs = std::filesystem;

PMP_File_Parse_Result read_json(
  const std::filesystem::path &input_path,
  const std::function<bool(size_t matrix_index)> &should_parse_matrix);

PMP_File_Parse_Result read_mathematica(
  const std::filesystem::path &input_path,
  const std::function<bool(size_t matrix_index)> &should_parse_matrix);

PMP_File_Parse_Result
read_xml(const std::filesystem::path &input_file,
         const std::function<bool(size_t matrix_index)> &should_parse_matrix);

PMP_File_Parse_Result PMP_File_Parse_Result::read(
  const fs::path &input_path,
  const std::function<bool(size_t matrix_index)> &should_parse_matrix)
{
  PMP_File_Parse_Result result;
  try
    {
      // Read

      if(input_path.extension() == ".json")
        {
          result = read_json(input_path, should_parse_matrix);
        }
      else if(input_path.extension() == ".m")
        {
          result = read_mathematica(input_path, should_parse_matrix);
        }
      else if(input_path.extension() == ".xml")
        {
          result = read_xml(input_path, should_parse_matrix);
        }
      else
        {
          El::RuntimeError("Expected .json, .m, or .xml extension.");
        }

      // Validate

      validate(result);
    }
  catch(std::exception &e)
    {
      El::RuntimeError("Error when parsing ", input_path, ": ", e.what());
    }
  return result;
}
void PMP_File_Parse_Result::validate(const PMP_File_Parse_Result &result)
{
  if(result.parsed_matrices.size() > result.num_matrices)
    El::LogicError("parsed_matrices.size()=", result.parsed_matrices.size(),
                   " should not exceed num_matrices=", result.num_matrices);

  for(auto &[index, matrix] : result.parsed_matrices)
    {
      if(index >= result.num_matrices)
        El::LogicError("index=", index, "should be less than",
                       result.num_matrices);
    }

  if(result.num_matrices == 0 && result.objective.empty()
     && result.normalization.empty())
    El::RuntimeError("Nothing was read from input file.");
}