#include "../parse_vector.hxx"
#include "../parse_number.hxx"
#include "../../../../Damped_Rational.hxx"

#include <boost/algorithm/string.hpp>

#include <algorithm>
#include <iterator>
#include <string>

std::vector<char>::const_iterator
parse_damped_rational(const std::vector<char>::const_iterator &begin,
                      const std::vector<char>::const_iterator &end,
                      Damped_Rational &damped_rational)
{
  const std::string damped_literal("DampedRational[");
  auto damped_start(
    std::search(begin, end, damped_literal.begin(), damped_literal.end()));
  if(damped_start == end)
    {
      throw std::runtime_error("Could not find '" + damped_literal + "'");
    }

  auto comma(std::find(damped_start, end, ','));
  if(comma == end)
    {
      throw std::runtime_error("Missing comma after DampedRational.constant");
    }
  auto constant_start(std::next(damped_start, damped_literal.size()));
  damped_rational.constant = Boost_Float(parse_number(constant_start, comma));

  auto start_poles(std::next(comma));
  auto end_poles(parse_vector(start_poles, end, damped_rational.poles));

  comma = std::find(end_poles, end, ',');
  if(comma == end)
    {
      throw std::runtime_error("Missing comma after DampedRational.poles");
    }

  auto start_base(std::next(comma));
  comma = std::find(start_base, end, ',');
  if(comma == end)
    {
      throw std::runtime_error("Missing comma after DampedRational.base");
    }
  damped_rational.base = Boost_Float(parse_number(start_base, comma));

  auto start_variable(std::next(comma));
  const auto close_bracket(std::find(start_variable, end, ']'));
  if(close_bracket == end)
    {
      throw std::runtime_error("Missing ']' at end of DampedRational");
    }
  damped_rational.variable
    = boost::trim_copy(std::string(start_variable, close_bracket));

  return std::next(close_bracket);
}
