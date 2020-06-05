#include "../Polynomial_State.hxx"

void Polynomial_State::json_key(const std::string &key)
{
  throw std::runtime_error("Invalid input file.  Unexpected key '" + key
                           + "' inside '" + name + "'.");
}
