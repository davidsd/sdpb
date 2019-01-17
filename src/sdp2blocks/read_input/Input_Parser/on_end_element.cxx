#include "../Input_Parser.hxx"

void Input_Parser::on_end_element(const std::string &element_name)
{
  if(inside_expression)
    {
      if(inside_sdp)
        {
          if(element_name != "Symbol")
            {
              if(objective_state.on_end_element(element_name))
                {
                  finished_objective = !objective_state.inside;
                }
              else if(normalization_state.on_end_element(element_name))
                {
                  finished_normalization = !normalization_state.inside;
                }
              else if(positive_matrix_with_prefactor_state.on_end_element(
                        element_name))
                {}
              else
                {
                  inside_sdp = (element_name != sdp_name);
                }
            }
        }
      else
        {
          inside_expression = (element_name != expression_name);
        }
    }
}
