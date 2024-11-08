#pragma once

#include "Damped_Rational.hxx"
#include "Polynomial_Power_Product.hxx"

#include <El.hpp>

#include <filesystem>
#include <optional>
#include <vector>

struct PVM_Info
{
  std::filesystem::path block_path;
  Damped_Rational prefactor;
  Damped_Rational reduced_prefactor;
  std::vector<El::BigFloat> sample_points;
  std::vector<El::BigFloat> sample_scalings;
  std::vector<El::BigFloat> reduced_sample_scalings;
  std::optional<std::vector<Polynomial_Power_Product>> preconditioning_vector;
};
