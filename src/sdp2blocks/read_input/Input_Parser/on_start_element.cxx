#include "../Input_Parser.hxx"

#include <stdexcept>

void Input_Parser::on_start_element(const std::string &element_name)
{
  if(inside_expression)
    {
      if(inside_sdp)
        {
          if(element_name != "Symbol"
             && !objective_state.on_start_element(element_name))
            {
              throw std::runtime_error(
                "Invalid input file.  Expected '" + objective_state.name
                + "' inside 'Expression.Function', but found '" + element_name
                + "'");
            }
        }
      else if(element_name == "Function")
        {
          inside_sdp = true;
        }
      else
        {
          throw std::runtime_error(
            "Invalid input file.  Expected 'Function' but found '"
            + element_name + "'");
        }
    }
  else if(element_name == expression_name)
    {
      inside_expression = true;
    }
  else
    {
      throw std::runtime_error(
        "Invalid input file.  Expected 'Expression' but found '" + element_name
        + "'");
    }
}
