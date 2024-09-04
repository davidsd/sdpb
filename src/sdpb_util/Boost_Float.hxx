#pragma once

#include <El.hpp>
#include <boost/multiprecision/mpfr.hpp>

using Boost_Float = boost::multiprecision::mpfr_float;

std::string to_string(const Boost_Float &boost_float);

Boost_Float to_Boost_Float(const El::BigFloat &alpha);

El::BigFloat to_BigFloat(const Boost_Float &value);

std::vector<El::BigFloat>
to_BigFloat_Vector(const std::vector<Boost_Float> &input);

std::vector<Boost_Float>
to_Boost_Float_Vector(const std::vector<El::BigFloat> &input);
