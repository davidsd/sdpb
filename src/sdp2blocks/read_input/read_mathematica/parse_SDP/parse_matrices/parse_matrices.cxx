#include "../../../Positive_Matrix_With_Prefactor_State.hxx"

#include <algorithm>
#include <iterator>
#include <string>

std::vector<char>::const_iterator
parse_matrix(const std::vector<char>::const_iterator &begin,
             const std::vector<char>::const_iterator &end,
             Positive_Matrix_With_Prefactor &matrix);

std::vector<char>::const_iterator
parse_matrices(const std::vector<char>::const_iterator &begin,
               const std::vector<char>::const_iterator &end,
               std::vector<Positive_Matrix_With_Prefactor> &matrices)
{
  const auto open_brace(std::find(begin, end, '{'));
  if(open_brace == end)
    {
      throw std::runtime_error("Could not find '{' to start array of "
                               "PositiveMatrixWithPrefactor");
    }

  auto delimiter(open_brace);
  const std::vector<char> delimiters({',', '}'});
  do
    {
      auto start_matrix(std::next(delimiter, 1));
      matrices.emplace_back();
      auto end_matrix(parse_matrix(start_matrix, end, matrices.back()));
      delimiter = std::find_first_of(end_matrix, end, delimiters.begin(),
                                     delimiters.end());
      if(delimiter == end)
        {
          throw std::runtime_error(
            "Missing '}' at end of array of PositiveMatrixWithPrefactor");
        }
    }
  while(*delimiter != '}');
}
