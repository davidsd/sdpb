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
        /* El::BuildString() is necessary in case of empty __VA_ARGS__ */    \
        RUNTIME_ERROR("Assertion '", #condition, "' failed: \n    ",          \
                      El::BuildString(__VA_ARGS__));                          \
  } while(false)
