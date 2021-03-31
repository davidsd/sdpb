#include "../Function_State.hxx"

void Function_State::json_start_object()
{
  if(inside)
    {
      throw std::runtime_error("Invalid input file.  Found an object inside '"
                               + name + "'.");
    }
  else
    {
      inside = true;
    }
}
