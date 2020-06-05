#include "../Polynomial_State.hxx"

void Polynomial_State::json_string(const std::string &s)
{
  vector_polynomial_term_state.json_string(s);
}
