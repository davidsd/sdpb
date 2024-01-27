#include "pmp_read/PMP_File_Parse_Result.hxx"
#include "sdpb_util/assert.hxx"

#include <boost/interprocess/file_mapping.hpp>
#include <boost/interprocess/mapped_region.hpp>

namespace fs = std::filesystem;

const char *
parse_SDP(const char *begin, const char *end,
          const std::function<bool(size_t matrix_index)> &should_parse_matrix,
          std::vector<El::BigFloat> &objectives,
          std::vector<El::BigFloat> &normalization, size_t &num_matrices,
          std::map<size_t, Polynomial_Vector_Matrix> &parsed_matrices);

PMP_File_Parse_Result read_mathematica(
  const std::filesystem::path &input_path,
  const std::function<bool(size_t matrix_index)> &should_parse_matrix)
{
  PMP_File_Parse_Result result;
  std::ifstream input_stream(input_path);
  if(!input_stream.good())
    {
      RUNTIME_ERROR("Unable to open input: ", input_path);
    }

  boost::interprocess::file_mapping mapped_file(
    input_path.c_str(), boost::interprocess::read_only);
  boost::interprocess::mapped_region mapped_region(
    mapped_file, boost::interprocess::read_only);

      const char *begin(
        static_cast<const char *>(mapped_region.get_address())),
        *end(begin + mapped_region.get_size());
      parse_SDP(begin, end, should_parse_matrix, result.objective,
                result.normalization, result.num_matrices,
                result.parsed_matrices);
  return result;
}
