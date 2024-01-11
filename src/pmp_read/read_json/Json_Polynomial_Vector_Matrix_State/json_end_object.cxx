#include "../Json_Polynomial_Vector_Matrix_State.hxx"
#include "sdpb_util/to_matrix.hxx"

namespace
{
  std::vector<Polynomial>
  to_poly_vector(const std::vector<std::vector<El::BigFloat>> &input_vectors)
  {
    std::vector<Polynomial> poly_vector;
    for(const auto &coeff_vector : input_vectors)
      {
        poly_vector.emplace_back().coefficients = coeff_vector;
      }
    return poly_vector;
  }
}

void Json_Polynomial_Vector_Matrix_State::json_end_object()
{
  if(parsing_damped_rational)
    {
      damped_rational_state.json_end_object();
      parsing_damped_rational = damped_rational_state.inside;
    }
  else if(parsing_polynomials)
    {
      throw std::runtime_error(
        "Invalid input file.  Found an object end inside '" + name + "."
        + polynomials_state.name + "'");
    }
  else
    {
      value = std::make_unique<Polynomial_Vector_Matrix>(
        to_matrix<std::vector<Polynomial>,
                  std::vector<std::vector<El::BigFloat>>>(
          polynomials_state.value, to_poly_vector),
        std::make_optional(damped_rational_state.value), std::nullopt,
        std::nullopt, std::nullopt);
      inside = false;
    }
}
