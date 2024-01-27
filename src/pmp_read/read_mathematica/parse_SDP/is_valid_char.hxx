#pragma once

#include <cctype>

inline bool is_valid_char(const char &c)
{
  return !std::isspace(c) && c != '\\';
}
