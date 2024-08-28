#pragma once

#include "Zero.hxx"
#include "pmp/Damped_Rational.hxx"
#include "pmp/Polynomial_Power_Product.hxx"

#include <optional>

void compute_lambda(const std::vector<El::BigFloat> &sample_points,
                    const std::vector<El::BigFloat> &sample_scalings,
                    const std::optional<std::vector<Polynomial_Power_Product>>
                      &preconditioning_vector,
                    const Damped_Rational &prefactor, const size_t &num_rows,
                    const El::Matrix<El::BigFloat> &x,
                    const std::vector<El::BigFloat> &zero_vector,
                    std::vector<Zero> &zeros, El::BigFloat &error);
