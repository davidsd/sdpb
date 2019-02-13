#pragma once

#include <El.hpp>

#include <vector>

El::BigFloat parse_number(const std::vector<char>::const_iterator &begin,
                          const std::vector<char>::const_iterator &end);
