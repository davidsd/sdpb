#include "Derivative_Term.hxx"

#include <set>

std::set<Derivative_Term> derivative(const Derivative_Term &term)
{
  // Handle the leading term of just multiplying by a single
  // derivative f'
  std::set<Derivative_Term> result;
  {
    Derivative_Term element(term);
    auto single_deriv_element(element.powers.find(1));
    if(single_deriv_element == element.powers.end())
      {
        element.powers.insert(std::make_pair(int64_t(1), int64_t(1)));
      }
    else
      {
        ++(single_deriv_element->second);
      }
    result.insert(element);
  }
  
  // Take derivative of all of the individual powers
  for(auto &power : term.powers)
    {
      Derivative_Term deriv(term);
      deriv.constant *= power.second;
      auto old_deriv(deriv.powers.find(power.first));
      --(old_deriv->second);
      if(old_deriv->second == 0)
        {
          deriv.powers.erase(old_deriv);
        }
      auto new_deriv(deriv.powers.find(power.first + 1));
      if(new_deriv == deriv.powers.end())
        {
          deriv.powers.insert(std::make_pair(power.first + 1, 1));
        }
      else
        {
          ++(new_deriv->second);
        }
      result.insert(deriv);
    }
  return result;
}
