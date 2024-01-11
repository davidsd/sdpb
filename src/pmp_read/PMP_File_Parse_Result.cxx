#include "PMP_File_Parse_Result.hxx"

#include "read_json/read_json.hxx"
#include "read_mathematica/read_mathematica.hxx"
#include "read_xml/read_xml.hxx"

namespace fs = std::filesystem;

PMP_File_Parse_Result::PMP_File_Parse_Result(
  const fs::path &input_path,
  const std::function<bool(size_t matrix_index)> &should_parse_matrix)
    : path(input_path), num_matrices(0)
{
  try
    {
      // Read

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
      else if(input_path.extension() == ".xml")
        {
          read_xml(input_path, should_parse_matrix, objective, num_matrices,
                   parsed_matrices);
        }
      else
        {
          El::RuntimeError("Expected .json, .m, or .xml extension.");
        }

      // Validate

      if(parsed_matrices.size() > num_matrices)
        El::LogicError("parsed_matrices.size()=", parsed_matrices.size(),
                       " should not exceed num_matrices=", num_matrices);

      for(auto &[index, matrix] : parsed_matrices)
        {
          if(index >= num_matrices)
            El::LogicError("index=", index, "should be less than",
                           num_matrices);
        }

      if(num_matrices == 0 && objective.empty() && normalization.empty())
        El::RuntimeError("Nothing was read from input file.");
    }
  catch(std::exception &e)
    {
      El::RuntimeError("Error when parsing '", input_path, "': ", e.what());
    }
}