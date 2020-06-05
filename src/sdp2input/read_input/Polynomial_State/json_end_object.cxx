#include "../Polynomial_State.hxx"

void Polynomial_State::json_end_object()
{
  throw std::runtime_error(
    "Invalid input file.  Unexpected object end inside '" + name + "'.");
}
