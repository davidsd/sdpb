#include "../parse_generic.hxx"
#include "sdpb_util/to_matrix.hxx"

#include <algorithm>
#include <iterator>
#include <string>

const char *parse_damped_rational(const char *begin, const char *end,
                                  Damped_Rational &damped_rational);

// TODO fill matrix
const char *
parse_matrix(const char *begin, const char *end,
             std::unique_ptr<Polynomial_Vector_Matrix> &matrix)
{
  const std::string matrix_literal("PositiveMatrixWithPrefactor[");
  auto matrix_start(
    std::search(begin, end, matrix_literal.begin(), matrix_literal.end()));
  if(matrix_start == end)
    {
      throw std::runtime_error("Could not find '" + matrix_literal + "'");
    }

  Damped_Rational prefactor;
  auto start_damped_rational(std::next(matrix_start, matrix_literal.size()));
  auto end_damped_rational(
    parse_damped_rational(start_damped_rational, end, prefactor));

  auto comma(std::find(end_damped_rational, end, ','));
  if(comma == end)
    {
      throw std::runtime_error("Missing comma after DampedRational");
    }

  std::vector<std::vector<std::vector<Polynomial>>> polynomials;
  const char *end_polynomials(
    parse_generic(std::next(comma), end, polynomials));

  const char *close_bracket(std::find(end_polynomials, end, ']'));
  if(close_bracket == end)
    {
      throw std::runtime_error("Missing ']' at end of SDP");
    }

  matrix = std::make_unique<Polynomial_Vector_Matrix>(
    to_matrix(polynomials), std::make_optional(prefactor), std::nullopt,
    std::nullopt, std::nullopt);
  return std::next(close_bracket);
}
