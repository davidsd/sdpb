#pragma once

#include <El.hpp>

#include <vector>

std::vector<char>::const_iterator
parse_vector(const std::vector<char>::const_iterator &begin,
             const std::vector<char>::const_iterator &end,
             std::vector<El::BigFloat> &result_vector);
