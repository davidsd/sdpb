#include "PVM_File_Parse_Result.hxx"

#include "read_xml/read_xml.hxx"

namespace fs = std::filesystem;

PVM_File_Parse_Result::PVM_File_Parse_Result(
  const fs::path &input_path,
  const std::function<bool(size_t matrix_index)> &should_parse_matrix)
    : path(input_path), num_matrices(0)
{
  if(input_path.extension() == ".xml")
    {
      read_xml(input_path, should_parse_matrix, objective, num_matrices, parsed_matrices);
    }
  else
    {
      El::RuntimeError("Cannot parse input file: ", input_path,
                       ". Expected .xml extension.");
    }

  for(auto &[index, matrix] : parsed_matrices)
    {
      if(index >= num_matrices)
        El::LogicError("index=", index, "should be less than", num_matrices,
                       ", input_path=", input_path);
    }
}