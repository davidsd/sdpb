#pragma once

#include <ostream>
#include <string>
#include <boost/optional.hpp>

struct Memory_Limit
{
  // Value passed to command line
  std::string str{};
  // Actual value in bytes
  size_t bytes{};
  Memory_Limit(const std::string &input, size_t bytes);

  Memory_Limit() = default;
  Memory_Limit(const Memory_Limit &other) = default;
  Memory_Limit(Memory_Limit &&other) noexcept = default;
  Memory_Limit &operator=(const Memory_Limit &other) = default;
  Memory_Limit &operator=(Memory_Limit &&other) noexcept = default;

  size_t limit_or_infinite() const;
  friend std::ostream &operator<<(std::ostream &os, const Memory_Limit &m);
};

// Converts string containing memory size in [kilo|mega|giga]bytes
// into number of bytes, e.g.:
// 100 or 100B -> 100
// 100K or 100KB -> 102400
// 100M or 100MB -> 104857600
// 100G or 100GB -> 107374182400
//
// Also supports percents of MemAvailable (measured at program start), e.g.
// 80% -> 0.8 * MemAvailable
//
// Based on example from
// https://theboostcpplibraries.com/boost.propertytree#ex.propertytree_03
struct String_To_Memory_Limit_Translator
{
  typedef std::string internal_type;
  typedef Memory_Limit external_type;

  const size_t mem_available_bytes;

  explicit String_To_Memory_Limit_Translator(size_t mem_available_bytes);

  Memory_Limit from_string(const std::string &s) const;

  boost::optional<Memory_Limit> get_value(const std::string &s) const;

  std::string put_value(const Memory_Limit &value);
};
