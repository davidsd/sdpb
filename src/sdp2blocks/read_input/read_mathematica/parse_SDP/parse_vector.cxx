#include "parse_number.hxx"

#include <iterator>

std::vector<char>::const_iterator
parse_vector(const std::vector<char>::const_iterator &begin,
             const std::vector<char>::const_iterator &end,
             std::vector<El::BigFloat> &result_vector)
{
  const auto open_brace(std::find(begin, end, '{'));
  if(open_brace == end)
    {
      throw std::runtime_error("Missing '{' at beginning of array of numbers");
    }
  auto start_element(std::next(open_brace, 1));

  const auto close_brace(std::find(start_element, end, '}'));
  if(close_brace == end)
    {
      throw std::runtime_error("Missing '}' at end of array of numbers");
    }

  auto comma(open_brace);
  comma = std::find(start_element, close_brace, ',');
  while(start_element < close_brace)
  {
    
    result_vector.emplace_back(parse_number(start_element, comma));
    start_element = std::next(comma, 1);
    comma = std::find(start_element, close_brace, ',');
  }

  return std::next(close_brace, 1);
}
