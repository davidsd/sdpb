#include <El.hpp>

#include <vector>
#include <cctype>
#include <string>

El::BigFloat parse_number(const std::vector<char>::const_iterator &begin,
                          const std::vector<char>::const_iterator &end)
{
  std::string cleaned_string;
  cleaned_string.reserve(end-begin);
  for(auto c(begin); c!=end; ++c)
    {
      if(!std::isspace(*c) && *c!='\\')
        {
          cleaned_string.push_back(*c);
        }
    }
  return El::BigFloat(cleaned_string);
}
