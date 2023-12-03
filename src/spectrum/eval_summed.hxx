#pragma once

#include "sdpb_util/Polynomial.hxx"

El::BigFloat
eval_summed(const std::vector<std::vector<Polynomial>> &summed_polynomials,
            const El::BigFloat &x);
