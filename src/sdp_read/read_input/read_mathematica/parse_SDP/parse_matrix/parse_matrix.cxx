#include "../parse_generic.hxx"
#include "sdp_read/Positive_Matrix_With_Prefactor.hxx"

#include <algorithm>
#include <iterator>
#include <string>

const char *parse_damped_rational(const char *begin, const char *end,
                                  Damped_Rational &damped_rational);

const char *parse_matrix(const char *begin, const char *end,
                         Positive_Matrix_With_Prefactor &matrix)
{
  const std::string matrix_literal("PositiveMatrixWithPrefactor[");
  auto matrix_start(
    std::search(begin, end, matrix_literal.begin(), matrix_literal.end()));
  if(matrix_start == end)
    {
      throw std::runtime_error("Could not find '" + matrix_literal + "'");
    }

  auto start_damped_rational(std::next(matrix_start, matrix_literal.size()));
  auto end_damped_rational(
    parse_damped_rational(start_damped_rational, end, matrix.damped_rational));

  auto comma(std::find(end_damped_rational, end, ','));
  if(comma == end)
    {
      throw std::runtime_error("Missing comma after DampedRational");
    }

  const char *end_polynomials(
    parse_generic(std::next(comma), end, matrix.polynomials));

  const char *close_bracket(std::find(end_polynomials, end, ']'));
  if(close_bracket == end)
    {
      throw std::runtime_error("Missing ']' at end of SDP");
    }
  return std::next(close_bracket);
}
