#include "../Function_State.hxx"

void Function_State::json_start_object()
{
  if(inside)
    {
      RUNTIME_ERROR("Invalid input file. Found an object inside '"
                               + name , "'.");
    }
  else
    {
      inside = true;
    }
}
