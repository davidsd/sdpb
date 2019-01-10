#include "../Input_Parser.hxx"

void Input_Parser::on_end_element(const std::string &element_name)
{
  if(inside_expression)
    {
      if(inside_sdp)
        {
          if(!finished_objective
             && objective_state.on_end_element(element_name))
            {
              finished_objective = true;
            }
          else if(!finished_normalization
                  && normalization_state.on_end_element(element_name))
            {
              finished_normalization = true;
            }
          else if(element_name == sdp_name)
            {
              inside_sdp = false;
            }
        }
      else if(element_name == expression_name)
        {
          inside_expression = false;
        }
    }
}
