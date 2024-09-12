#pragma once

#include "pmp/Damped_Rational.hxx"
#include "Polynomial_Matrix_Program.hxx"

#include <filesystem>
#include <utility>
#include <vector>

struct PVM_Info
{
  int block_index = -1;
  std::filesystem::path block_path;
  size_t dim = 0;
  Damped_Rational prefactor;
  Damped_Rational reduced_prefactor;
  std::vector<El::BigFloat> sample_points;
  std::vector<El::BigFloat> sample_scalings;
  std::vector<El::BigFloat> reduced_sample_scalings;

  PVM_Info() = default;

  PVM_Info(const Polynomial_Vector_Matrix &pvm, int block_index,
           std::filesystem::path block_path)
      : block_index(block_index),
        block_path(std::move(block_path)),
        dim(pvm.polynomials.Height()),
        prefactor(pvm.prefactor),
        reduced_prefactor(pvm.reduced_prefactor),
        sample_points(pvm.sample_points),
        sample_scalings(pvm.sample_scalings),
        reduced_sample_scalings(pvm.reduced_sample_scalings)
  {}

  void validate() const
  {
    ASSERT(block_index >= 0, DEBUG_STRING(block_index));
    ASSERT(!block_path.empty());
    ASSERT(dim > 0);
    const auto num_points = sample_points.size();
    ASSERT_EQUAL(num_points, sample_scalings.size());
    ASSERT_EQUAL(num_points, reduced_sample_scalings.size());
  }

  void validate(const size_t num_blocks) const
  {
    validate();
    ASSERT(block_index < num_blocks, DEBUG_STRING(block_index),
           DEBUG_STRING(num_blocks));
  }
};

struct PMP_Info
{
  // Total number of blocks
  size_t num_blocks;
  // Blocks stored on the current rank
  std::vector<PVM_Info> blocks;

  explicit
  PMP_Info(const size_t num_blocks, const std::vector<PVM_Info> &blocks)
      : num_blocks(num_blocks), blocks(blocks)
  {}

  explicit PMP_Info(const Polynomial_Matrix_Program &pmp)
      : num_blocks(pmp.num_matrices)
  {
    const auto num_local_blocks = pmp.matrix_index_local_to_global.size();
    blocks.reserve(num_local_blocks);
    for(size_t i = 0; i < num_local_blocks; ++i)
      {
        blocks.emplace_back(pmp.matrices.at(i),
                            pmp.matrix_index_local_to_global.at(i),
                            pmp.block_paths.at(i));
      }
  }
};
