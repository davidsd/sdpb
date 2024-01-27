#include "is_valid_char.hxx"

#include <algorithm>
#include <vector>
#include <string>

std::string parse_number(const char *begin, const char *end)
{
  std::string cleaned_string;
  cleaned_string.reserve(end - begin);
  auto c(begin);
  for(; c != end && *c != '`'; ++c)
    {
      if(is_valid_char(*c))
        {
          cleaned_string.push_back(*c);
        }
    }
  auto carat(std::find(c, end, '^'));
  if(carat != end)
    {
      cleaned_string.push_back('e');
      for(auto c(std::next(carat)); c != end; ++c)
        {
          if(is_valid_char(*c))
            {
              cleaned_string.push_back(*c);
            }
        }
    }
  return cleaned_string;
}
