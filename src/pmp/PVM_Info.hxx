#pragma once

#include "Damped_Rational.hxx"

#include <El.hpp>

#include <filesystem>
#include <vector>

struct PVM_Info
{
  std::filesystem::path block_path;
  Damped_Rational prefactor;
  Damped_Rational reduced_prefactor;
  std::vector<El::BigFloat> sample_points;
  std::vector<El::BigFloat> sample_scalings;
  std::vector<El::BigFloat> reduced_sample_scalings;
};
