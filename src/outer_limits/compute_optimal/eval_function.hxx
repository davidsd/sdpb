#pragma once

#include <El.hpp>

#include <map>

El::BigFloat
eval_function(const El::BigFloat &infinity,
              const std::map<El::BigFloat, El::BigFloat> &f, const El::BigFloat &x);
