#include "../Positive_Matrix_With_Prefactor_State.hxx"

bool Positive_Matrix_With_Prefactor_State::xml_on_end_element(
  const std::string &element_name)
{
  bool result(inside);
  if(inside)
    {
      if(damped_rational_state.xml_on_end_element(element_name))
        {
          finished_damped_rational = !damped_rational_state.inside;
          if(finished_damped_rational)
            {
              using namespace std;
              swap(damped_rational_state.value, value.damped_rational);
            }
        }
      else if(polynomials_state.xml_on_end_element(element_name))
        {
          if(!polynomials_state.inside)
            {
              value.polynomials.reserve(polynomials_state.value.size());
              for(auto &polynomial_vector_vector : polynomials_state.value)
                {
                  value.polynomials.emplace_back();
                  auto &vector_vector_poly(value.polynomials.back());
                  vector_vector_poly.reserve(polynomial_vector_vector.size());
                  for(auto &polynomial_vector : polynomial_vector_vector)
                    {
                      vector_vector_poly.emplace_back();
                      auto &vector_poly(vector_vector_poly.back());
                      vector_poly.reserve(polynomial_vector.size());
                      for(auto &polynomial : polynomial_vector)
                        {
                          vector_poly.emplace_back();
                          auto &poly(vector_poly.back());
                          poly.coefficients.resize(polynomial.size(), 0);
                          for(auto &term : polynomial)
                            {
                              if(poly.coefficients.size() < term.first + 1)
                                {
                                  poly.coefficients.resize(term.first + 1, 0);
                                }
                              poly.coefficients[term.first] = term.second;
                            }
                        }
                    }
                }
            }
        }
      else
        {
          inside = (element_name != name);
        }
    }
  return result;
}
