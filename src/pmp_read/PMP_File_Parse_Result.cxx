#include "PMP_File_Parse_Result.hxx"

#include "sdpb_util/assert.hxx"

#include <boost/algorithm/string/predicate.hpp>

namespace fs = std::filesystem;

PMP_File_Parse_Result
read_json(const std::filesystem::path &input_path, bool should_parse_objective,
          bool should_parse_normalization,
          const std::function<bool(size_t matrix_index)> &should_parse_matrix);

PMP_File_Parse_Result read_mathematica(
  const std::filesystem::path &input_path,
  const std::function<bool(size_t matrix_index)> &should_parse_matrix);

PMP_File_Parse_Result
read_xml(const std::filesystem::path &input_file,
         const std::function<bool(size_t matrix_index)> &should_parse_matrix);

PMP_File_Parse_Result
read_sdpa(const std::filesystem::path &input_file,
          const std::function<bool(size_t matrix_index)> &should_parse_matrix);

PMP_File_Parse_Result PMP_File_Parse_Result::read(
  const fs::path &input_path, bool should_parse_objective,
  bool should_parse_normalization,
  const std::function<bool(size_t matrix_index)> &should_parse_matrix)
{
  PMP_File_Parse_Result result;
  try
    {
      // Read

      if(input_path.extension() == ".json")
        {
          result = read_json(input_path, should_parse_objective,
                             should_parse_normalization, should_parse_matrix);
        }
      // TODO: use should_parse_objective and should_parse_normalization
      // also in read_mathematica() and read_xml()
      else if(input_path.extension() == ".m")
        {
          result = read_mathematica(input_path, should_parse_matrix);
        }
      else if(input_path.extension() == ".xml")
        {
          result = read_xml(input_path, should_parse_matrix);
        }
      // .dat: SDPA format
      // .dat-s: sparse SDPA format, see e.g. https://github.com/vsdp/SDPLIB
      // TODO support also CBF (Conic Becnhmark Format) https://cblib.zib.de/
      else if(input_path.extension() == ".dat"
              || input_path.extension() == ".dat-s"
              || boost::algorithm::ends_with(input_path.string(), ".dat.gz")
              || boost::algorithm::ends_with(input_path.string(), ".dat-s.gz"))
        {
          result = read_sdpa(input_path, should_parse_matrix);
        }
      else
        {
          RUNTIME_ERROR("Expected .json, .m, .xml, .dat or .dat-s extension.",
                        DEBUG_STRING(input_path));
        }

      // Validate

      validate(result);
    }
  catch(std::exception &e)
    {
      RUNTIME_ERROR("Error when parsing ", input_path, ": ", e.what());
    }
  return result;
}
void PMP_File_Parse_Result::validate(const PMP_File_Parse_Result &result)
{
  ASSERT(result.parsed_matrices.size() <= result.num_matrices,
         DEBUG_STRING(result.parsed_matrices.size()),
         DEBUG_STRING(result.num_matrices));

  for(auto &[index, matrix] : result.parsed_matrices)
    {
      ASSERT(index < result.num_matrices, DEBUG_STRING(index),
             DEBUG_STRING(result.num_matrices));
    }

  ASSERT(result.num_matrices > 0 || result.objective.has_value()
           || result.normalization.has_value(),
         "Nothing was read from input file.");
}