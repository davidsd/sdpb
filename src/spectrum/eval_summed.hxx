#pragma once

#include "pmp/Polynomial.hxx"

El::BigFloat
eval_summed(const std::vector<Polynomial_Vector> &summed_polynomials,
            const El::BigFloat &x);
