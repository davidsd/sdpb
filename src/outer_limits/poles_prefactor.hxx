#pragma once

#include "sdpb_util/Boost_Float.hxx"

#include <El.hpp>

#include <vector>

El::BigFloat
poles_prefactor(const std::vector<Boost_Float> &poles, const El::BigFloat &x);
