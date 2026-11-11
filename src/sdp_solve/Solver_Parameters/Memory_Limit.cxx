#include "Memory_Limit.hxx"

#include "sdpb_util/assert.hxx"
#include "sdpb_util/ostream/pretty_print_bytes.hxx"

#include <sstream>

Memory_Limit::Memory_Limit(const std::string &input, const size_t bytes)
    : str(input), bytes(bytes)
{}
size_t Memory_Limit::limit_or_infinite() const
{
  if(bytes == 0)
    return std::numeric_limits<size_t>::max();
  return bytes;
}
String_To_Memory_Limit_Translator::String_To_Memory_Limit_Translator(
  const size_t mem_available_bytes)
    : mem_available_bytes(mem_available_bytes)
{}
Memory_Limit
String_To_Memory_Limit_Translator::from_string(const std::string &s) const
{
  std::istringstream iss(s);
  double number;
  std::string suffix;

  iss >> number >> suffix;

  double factor;
  if(suffix.empty() || suffix == "B")
    factor = 1;
  else if(suffix == "K" || suffix == "KB" || suffix == "KiB")
    factor = 1024;
  else if(suffix == "M" || suffix == "MB" || suffix == "MiB")
    factor = 1024 * 1024;
  else if(suffix == "G" || suffix == "GB" || suffix == "GiB")
    factor = 1024 * 1024 * 1024;
  else if(suffix == "%")
    {
      ASSERT(number == 0 || mem_available_bytes != 0,
             "MemAvailable=0, cannot compute memory limit = ", s);
      factor = 0.01 * mem_available_bytes;
    }
  else
    RUNTIME_ERROR("Cannot parse memory size: \"", s, "\"");

  return Memory_Limit(s, number * factor);
}
boost::optional<Memory_Limit>
String_To_Memory_Limit_Translator::get_value(const std::string &s) const
{
  return from_string(s);
}
std::string
String_To_Memory_Limit_Translator::put_value(const Memory_Limit &value)
{
  return value.str;
}
std::ostream &operator<<(std::ostream &os, const Memory_Limit &m)
{
  return os << m.str << " (" << pretty_print_bytes(m.bytes) << ")";
}