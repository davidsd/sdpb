#include <algorithm>
#include <vector>
#include <cctype>
#include <string>

namespace
{
  inline bool is_valid_char(const char &c)
  {
    return !std::isspace(c) && c != '\\';
  }
}

std::string parse_number(const std::vector<char>::const_iterator &begin,
                         const std::vector<char>::const_iterator &end)
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
      for(auto c(std::next(carat, 1)); c != end; ++c)
        {
          if(is_valid_char(*c))
            {
              cleaned_string.push_back(*c);
            }
        }
    }
  return cleaned_string;
}
