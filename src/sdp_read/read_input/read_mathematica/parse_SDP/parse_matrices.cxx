#include "parse_vector.hxx"
#include "parse_generic.hxx"
#include "sdp_read/Positive_Matrix_With_Prefactor.hxx"

#include <algorithm>
#include <iterator>
#include <string>

const char *
parse_matrices(const char *begin, const char *end, const int &rank,
               const int &num_procs, const size_t &num_matrices,
               std::vector<Positive_Matrix_With_Prefactor> &matrices)
{
  const auto open_brace(std::find(begin, end, '{'));
  if(open_brace == end)
    {
      throw std::runtime_error("Could not find '{' to start array");
    }

  auto delimiter(open_brace);
  const std::vector<char> delimiters({',', '}'});
  int matrix_index(num_matrices);
  do
    {
      auto start_matrix(std::next(delimiter));
      Positive_Matrix_With_Prefactor matrix;
      auto end_matrix(parse_generic(start_matrix, end, matrix));
      matrices.emplace_back();
      if(matrix_index % num_procs == rank)
        {
          swap(matrices.back(), matrix);
        }
      ++matrix_index;

      delimiter = std::find_first_of(end_matrix, end, delimiters.begin(),
                                     delimiters.end());
      if(delimiter == end)
        {
          throw std::runtime_error(
            "Missing '}' at end of array of PositiveMatrixWithPrefactor");
        }
    }
  while(*delimiter != '}');
  return std::next(delimiter);
}
