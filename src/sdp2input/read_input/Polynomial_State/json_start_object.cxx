#include "../Polynomial_State.hxx"

void Polynomial_State::json_start_object()
{
  throw std::runtime_error("Invalid input file.  Unexpected object inside '"
                           + name + "'.");
}
