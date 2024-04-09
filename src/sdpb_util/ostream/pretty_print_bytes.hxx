#pragma once

#include <El.hpp>
#include <iomanip>

// Pretty pring memory size given in bytes.
// Examples:
// also_print_exact=false:
// - 100 B
// - 3.1 KB
// - 1250.0 GB
// also_print_exact=true:
// - 100 B
// - 3.1 KB (3180 bytes)
inline std::string
pretty_print_bytes(const size_t bytes, const bool also_print_exact = false)
{
  if(bytes < 1024)
    return El::BuildString(bytes, " B");
  double value = bytes;

  std::string unit;
  for(const auto curr_unit : {"KB", "MB", "GB"})
    {
      value /= 1024;
      unit = curr_unit;
      if(value < 1024)
        break;
    }

  std::ostringstream os;
  El::BuildStream(os, std::fixed, std::setprecision(1), value, " ", unit);
  if(also_print_exact)
    El::BuildStream(os, " (", bytes, " bytes)");
  return os.str();
}
