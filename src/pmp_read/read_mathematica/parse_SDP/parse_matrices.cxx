#include "parse_generic.hxx"
#include "pmp/Polynomial_Vector_Matrix.hxx"
#include "sdpb_util/to_matrix.hxx"

#include <algorithm>
#include <iterator>

// should_parse_matrix() accepts matrix index in file (starting from 0)
// and tells if this rank should parse the matrix of skip it.
const char *parse_matrices(
  const char *begin, const char *end,
  const std::function<bool(size_t matrix_index)> &should_parse_matrix,
  size_t &num_matrices,
  std::map<size_t, Polynomial_Vector_Matrix> &parsed_matrices)
{
  const auto open_brace(std::find(begin, end, '{'));
  if(open_brace == end)
    {
      RUNTIME_ERROR("Could not find '{' to start array");
    }

  auto delimiter(open_brace);
  const std::vector<char> delimiters({',', '}'});
  int matrix_index = 0;
  parsed_matrices.clear();
  do
    {
      auto start_matrix(std::next(delimiter));
      // TODO fast-forward through matrices that we don't need
      std::unique_ptr<Polynomial_Vector_Matrix> matrix;
      auto end_matrix(parse_matrix(start_matrix, end, matrix));
      if(should_parse_matrix(matrix_index))
        {
          parsed_matrices.emplace(matrix_index, std::move(*matrix));
        }
      ++matrix_index;

      delimiter = std::find_first_of(end_matrix, end, delimiters.begin(),
                                     delimiters.end());
      if(delimiter == end)
        {
          RUNTIME_ERROR(
            "Missing '}' at end of array of PositiveMatrixWithPrefactor");
        }
  } while(*delimiter != '}');
  num_matrices = matrix_index;
  return std::next(delimiter);
}
