#include "../XML_Parser.hxx"

#include <stdexcept>

void XML_Parser::on_start_element(const std::string &element_name)
{
  if(inside_expression)
    {
      if(inside_sdp)
        {
          if(element_name != "Symbol"
             && (finished_objective
                 || !objective_state.on_start_element(element_name))
             && (finished_normalization
                 || !normalization_state.on_start_element(element_name))
             && !positive_matrices_with_prefactor_state.on_start_element(
                  element_name))
            {
              throw std::runtime_error(
                "Invalid input file.  Expected 'Function' inside "
                "'Expression.Function', but found '"
                + element_name + "'");
            }
        }
      else
        {
          inside_sdp = (element_name == sdp_name);
          if(!inside_sdp)
            {
              throw std::runtime_error(
                "Invalid input file.  Expected 'Function' inside "
                "'Expression', but found '"
                + element_name + "'");
            }
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
