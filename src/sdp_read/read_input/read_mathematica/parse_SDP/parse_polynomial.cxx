#include "is_valid_char.hxx"
#include "sdp_read/Positive_Matrix_With_Prefactor.hxx"

#include <algorithm>
#include <iterator>
#include <string>

using namespace std::literals;

inline void check_iterator(const char character, const char *begin,
                           const char *iterator, const char *end)
{
  if(iterator == end)
    {
      throw std::runtime_error("Invalid polynomial string after '"s + character
                               + "': " + std::string(begin, end));
    }
}

const char *
parse_polynomial(const char *begin, const char *end, Polynomial &polynomial)
{
  const std::vector<char> delimiters({',', '}'});
  const auto delimiter(
    std::find_first_of(begin, end, delimiters.begin(), delimiters.end()));
  if(delimiter == end)
    {
      throw std::runtime_error("Missing '}' at end of array of polynomials");
    }

  std::string mantissa;
  for(auto c(begin); c < delimiter; ++c)
    {
      if(*c == '`')
        {
          do
            {
              ++c;
          } while(c != delimiter
                  && (std::isdigit(*c) || *c == '.' || !is_valid_char(*c)
                      || *c == '`'));
        }
      // Coefficients can be a plain 'x' instead of '1.0*x'.
      if(*c == '*' || *c == 'x')
        {
          std::string exponent;
          if(*c == '*')
            {
              do
                {
                  ++c;
                  check_iterator('*', begin, c, delimiter);
              } while(!is_valid_char(*c));

              if(*c == '^')
                {
                  exponent = "E";
                  ++c;
                  check_iterator('^', begin, c, delimiter);
                  while(c != delimiter
                        && ((exponent.size() == 1 && (*c == '-' || *c == '+'))
                            || std::isdigit(*c) || !is_valid_char(*c)))
                    {
                      if(is_valid_char(*c))
                        {
                          exponent.push_back(*c);
                        }
                      ++c;
                    }
                  while(c != delimiter && (!is_valid_char(*c) || *c == '*'))
                    {
                      ++c;
                    }
                }
            }
          else
            {
              mantissa = "1";
            }
          size_t degree(0);
          // Hard code the polynomial to be in 'x' since that is what
          // SDPB.m uses.
          if(*c == 'x')
            {
              ++c;
              while(!is_valid_char(*c))
                {
                  ++c;
                  check_iterator('x', begin, c, delimiter);
                }
              if(*c != '^')
                {
                  degree = 1;
                }
              else
                {
                  ++c;
                  std::string degree_string;
                  while((degree_string.empty() && (*c == '-' || *c == '+'))
                        || std::isdigit(*c) || !is_valid_char(*c))
                    {
                      if(is_valid_char(*c) && *c != '+')
                        {
                          degree_string.push_back(*c);
                        }
                      ++c;
                    }
                  degree = std::stoull(degree_string);
                }
            }

          if(polynomial.coefficients.size() < degree + 1)
            {
              polynomial.coefficients.resize(degree + 1);
            }
          polynomial.coefficients.at(degree)
            = El::BigFloat(mantissa + exponent);
          mantissa.clear();
        }
      else if(!mantissa.empty() && (*c == '-' || *c == '+' || c == delimiter))
        {
          if(polynomial.coefficients.size() < 1)
            {
              polynomial.coefficients.resize(1);
            }
          polynomial.coefficients.at(0) = El::BigFloat(mantissa);
          mantissa.clear();
        }
      if(c != delimiter && is_valid_char(*c) && *c != '+')
        {
          mantissa.push_back(*c);
        }
    }
  if(!mantissa.empty())
    {
      if(polynomial.coefficients.size() < 1)
        {
          polynomial.coefficients.resize(1);
        }
      polynomial.coefficients.at(0) = El::BigFloat(mantissa);
    }
  return delimiter;
}
