#pragma once

#include <El.hpp>

#include <exception>

#define THROW(exception_type, ...)                                            \
  do                                                                          \
    {                                                                         \
      throw exception_type(El::BuildString("in ", __FUNCTION__, "() at ",     \
                                           __FILE__, ":", __LINE__, ": \n  ", \
                                           __VA_ARGS__));                     \
  } while(false)

#define RUNTIME_ERROR(...) THROW(std::runtime_error, __VA_ARGS__)

#define LOGIC_ERROR(...) THROW(std::logic_error, __VA_ARGS__)

#define ASSERT(condition, ...)                                                \
  do                                                                          \
    {                                                                         \
      if(!(condition))                                                        \
        /* El::BuildString() is necessary in case of empty __VA_ARGS__ */     \
        RUNTIME_ERROR("Assertion '", #condition, "' failed:\n    ",           \
                      El::BuildString(__VA_ARGS__));                          \
  } while(false)

// Example:
// int x = 1;
// auto s = DEBUG_STRING(x+x) // s = "x+x=='2' "
// Trailing whitespace added for convenience
// in cases like ASSERT(false, DEBUG_STRING(a), DEBUG_STRING(b))
#define DEBUG_STRING(expr) El::BuildString(#expr, "='", expr, "' ")

// Example:
// int x = 1;
// int y = 2;
// ASSERT_EQUAL(x,y) // Assertion 'x == y' failed: x='1' x='2'
#define ASSERT_EQUAL(a, b, ...)                                               \
  ASSERT(a == b, DEBUG_STRING(a), "\n    ", DEBUG_STRING(b), "\n    ",        \
         El::BuildString(__VA_ARGS__))

#define PRINT_WARNING(...)                                                    \
  std::cerr << El::BuildString("Warning: ", __VA_ARGS__, "\n")
