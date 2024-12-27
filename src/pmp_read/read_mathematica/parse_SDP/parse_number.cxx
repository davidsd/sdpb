#include "is_valid_char.hxx"
#include "sdpb_util/assert.hxx"

#include <string>

const char *next_valid_char(const char *c, const char *end)
{
  do
    {
      ++c;
      if(is_valid_char(*c))
        break;
  } while(c != end);
  return c;
}

// "-1.234`8.45 *^ -10" -> "-1.234e-10"
std::string parse_number(const char *begin, const char *end)
{
  std::string cleaned_string;
  cleaned_string.reserve(end - begin);

  auto c = begin;

  // Go until precision mark "`" or exponent "*^"
  for(; c != end; c = next_valid_char(c, end))
    {
      if(*c == '`' || *c == '*')
        break;
      cleaned_string.push_back(*c);
    }

  // Skip precision mark:
  if(*c == '`')
    c = std::find(c, end, '*');
  if(c == end)
    return cleaned_string;

  // Replace "*^" with "e"
  c = next_valid_char(c, end);
  ASSERT_EQUAL(*c, '^');
  c = next_valid_char(c, end);
  cleaned_string.push_back('e');

  // Add exponent part
  for(; c != end; c = next_valid_char(c, end))
    {
      cleaned_string.push_back(*c);
    }

  return cleaned_string;
}
