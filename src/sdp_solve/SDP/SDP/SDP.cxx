#include "assign_bilinear_bases_dist.hxx"
#include "set_bases_blocks.hxx"
#include "sdp_solve/SDP.hxx"

#include "sdpb_util/assert.hxx"
#include "sdpb_util/copy_matrix.hxx"

#include <filesystem>

#include <iterator>

namespace fs = std::filesystem;

void read_block_data(const fs::path &sdp_path, const El::Grid &grid,
                     const Block_Info &block_info, SDP &sdp, Timers &timers);
void read_objectives(const fs::path &sdp_path, const El::Grid &grid,
                     El::BigFloat &objective_const,
                     El::DistMatrix<El::BigFloat> &dual_objective_b,
                     Timers &timers);

SDP::SDP(const fs::path &sdp_path, const Block_Info &block_info,
         const El::Grid &grid, Timers &timers)
{
  Scoped_Timer timer(timers, "SDP_ctor");
  read_objectives(sdp_path, grid, objective_const, dual_objective_b, timers);
  read_block_data(sdp_path, grid, block_info, *this, timers);

  Scoped_Timer validate_timer(timers, "validate");
  validate(block_info);
}

SDP::SDP(
  const El::BigFloat &objective_const_input,
  const std::vector<std::vector<El::BigFloat>> &primal_objective_c_input,
  const std::vector<El::Matrix<El::BigFloat>> &free_var_input,
  const El::DistMatrix<El::BigFloat, El::STAR, El::STAR> &yp_to_y_star,
  const El::DistMatrix<El::BigFloat, El::STAR, El::STAR> &dual_objective_b_star,
  const El::BigFloat &primal_c_scale, const Block_Info &block_info,
  const El::Grid &grid)
    : dual_objective_b(dual_objective_b_star.Height(), 1, grid),
      objective_const(objective_const_input)
{
  std::vector<size_t> block_offsets(primal_objective_c_input.size() + 1, 0);
  for(size_t p(0); p < primal_objective_c_input.size(); ++p)
    {
      block_offsets[p + 1]
        = block_offsets[p] + primal_objective_c_input[p].size();
    }

  auto &block_indices(block_info.block_indices);
  std::vector<El::Matrix<El::BigFloat>> bilinear_bases_local(
    2 * block_indices.size());
  for(size_t block(0); block != block_indices.size(); ++block)
    {
      bilinear_bases_local[2 * block].Resize(1, 1);
      bilinear_bases_local[2 * block + 1].Resize(0, 1);
      bilinear_bases_local[2 * block](0, 0) = 1;
    }
  assign_bilinear_bases_dist(bilinear_bases_local, grid, bilinear_bases);
  set_bases_blocks(block_info, bilinear_bases_local, bases_blocks, grid);

  for(size_t block(0); block != block_indices.size(); ++block)
    {
      primal_objective_c.blocks.emplace_back(
        primal_objective_c_input.at(block_indices.at(block)).size(), 1, grid);
      // TODO: copied from sdp_solve/SDP/SDP/read_primal_objective_c.cxx
      auto &primal_block(primal_objective_c.blocks.back());
      const size_t local_height(primal_block.LocalHeight());
      if(primal_block.GlobalCol(0) == 0)
        {
          for(size_t row = 0; row < local_height; ++row)
            {
              const size_t global_row(primal_block.GlobalRow(row));
              primal_block.SetLocal(
                row, 0,
                primal_c_scale
                  * primal_objective_c_input.at(block_indices.at(block))
                      .at(global_row));
            }
        }
    }

  // Setup B
  const int64_t B_Height(block_offsets.back()),
    B_Width(dual_objective_b_star.Height());
  El::Grid global_grid;
  El::DistMatrix<El::BigFloat> B(B_Height, B_Width, global_grid);

  for(int64_t row(0); row != B.LocalHeight(); ++row)
    {
      const size_t global_row(B.GlobalRow(row));
      const auto upper_iterator(std::upper_bound(
        block_offsets.begin(), block_offsets.end(), global_row));
      const size_t block_index(
        std::distance(block_offsets.begin(), upper_iterator) - 1);
      const size_t block_row(global_row - block_offsets[block_index]);
      for(int64_t column(0); column != B.LocalWidth(); ++column)
        {
          const size_t global_column(B.GlobalCol(column));
          B.SetLocal(
            row, column,
            primal_c_scale
              * free_var_input.at(block_index)(block_row, global_column));
        }
    }

  // Transform B into the yp frame.
  El::DistMatrix<El::BigFloat> U(B.Height(), B.Width(), B.Grid());
  El::Zero(U);
  El::DistMatrix<El::BigFloat> yp_to_y(yp_to_y_star.Height(),
                                       yp_to_y_star.Width(), B.Grid());

  copy_matrix(yp_to_y_star, yp_to_y);
  El::Gemm(El::Orientation::NORMAL, El::Orientation::NORMAL, El::BigFloat(1.0),
           B, yp_to_y, El::BigFloat(0.0), U);

  El::DistMatrix<El::BigFloat, El::STAR, El::STAR> U_star(U);

  free_var_matrix.blocks.reserve(block_indices.size());
  for(size_t block(0); block != block_indices.size(); ++block)
    {
      const size_t block_index(block_indices.at(block));
      free_var_matrix.blocks.emplace_back(
        free_var_input.at(block_index).Height(), U_star.Width(), grid);
      auto &free_var_block(free_var_matrix.blocks.back());
      for(int64_t row(0); row != free_var_block.Height(); ++row)
        {
          const int64_t global_row(block_offsets.at(block_index) + row);
          for(int64_t column(0); column != free_var_block.Width(); ++column)
            {
              if(free_var_block.IsLocal(row, column))
                {
                  free_var_block.SetLocal(free_var_block.LocalRow(row),
                                          free_var_block.LocalCol(column),
                                          U_star.GetLocal(global_row, column));
                }
            }
        }
    }

  // Copy over dual_objective_b
  copy_matrix(dual_objective_b_star, dual_objective_b);

  validate(block_info);
}

void SDP::validate(const Block_Info &block_info) const noexcept(false)
{
  const size_t num_blocks = block_info.block_indices.size();

  const auto error_prefix = "Invalid SDP: ";

  {
    // Check array sizes
    ASSERT_EQUAL(primal_objective_c.blocks.size(), num_blocks, error_prefix);
    ASSERT_EQUAL(free_var_matrix.blocks.size(), num_blocks, error_prefix);
    ASSERT_EQUAL(bilinear_bases.size(), 2 * num_blocks, error_prefix);
    ASSERT_EQUAL(bases_blocks.size(), 2 * num_blocks);
  }

  // Check dual_objective_b
  const auto &b = dual_objective_b;
  {
    ASSERT(El::mpi::Congruent(b.DistComm(), block_info.mpi_comm.value),
           error_prefix, " wrong MPI communicator for dual_objective_b");
    ASSERT_EQUAL(b.Width(), 1, error_prefix);
  }

  for(size_t index = 0; index < num_blocks; ++index)
    {
      const auto block_index = block_info.block_indices.at(index);

      const auto error_prefix_index
        = El::BuildString(error_prefix, "block_index = ", block_index, ": ");

      {
        // Check vector c and matrix B
        const auto &c = primal_objective_c.blocks.at(index);
        const auto &B = free_var_matrix.blocks.at(index);

        ASSERT(El::mpi::Congruent(c.DistComm(), block_info.mpi_comm.value),
               error_prefix_index,
               " wrong MPI communicator for primal_objective_c");
        ASSERT(El::mpi::Congruent(B.DistComm(), block_info.mpi_comm.value),
               error_prefix_index,
               " wrong MPI communicator for free_var_matrix");

        const auto P = block_info.get_schur_block_size(block_index);
        const auto N = b.Height();

        ASSERT_EQUAL(c.Height(), P, error_prefix_index);
        ASSERT_EQUAL(c.Width(), 1, error_prefix_index);
        ASSERT_EQUAL(B.Height(), P, error_prefix_index);
        ASSERT_EQUAL(B.Width(), N, error_prefix_index);
      }

      for(const size_t parity : {0, 1})
        {
          const auto error_prefix_index_parity
            = El::BuildString(error_prefix_index, "parity = ", parity, ": ");

          {
            // Check bilinear_bases matrix size
            // Should be consistent with Dual_Constraint_Group::Dual_Constraint_Group()
            auto &bilinear_bases_matrix
              = bilinear_bases.at(2 * index + parity);

            ASSERT(El::mpi::Congruent(bilinear_bases_matrix.DistComm(),
                                      block_info.mpi_comm.value),
                   error_prefix_index_parity,
                   " wrong MPI communicator for bilinear_bases");

            const size_t height
              = block_info.get_bilinear_bases_height(block_index, parity);
            const size_t width
              = block_info.get_bilinear_bases_width(block_index, parity);

            ASSERT_EQUAL(bilinear_bases_matrix.Height(), height,
                         error_prefix_index_parity);
            ASSERT_EQUAL(bilinear_bases_matrix.Width(), width,
                         error_prefix_index_parity);
          }

          {
            // Check bases_block size
            const auto &bases_block = bases_blocks.at(2 * index + parity);

            ASSERT(El::mpi::Congruent(bases_block.DistComm(),
                                      block_info.mpi_comm.value),
                   error_prefix_index_parity,
                   " wrong MPI communicator for bases_blocks");

            const auto height
              = block_info.get_psd_matrix_block_size(block_index, parity);
            const auto width = block_info.get_bilinear_pairing_block_size(
              block_index, parity);

            ASSERT_EQUAL(bases_block.Height(), height,
                         error_prefix_index_parity);
            ASSERT_EQUAL(bases_block.Width(), width,
                         error_prefix_index_parity);
          }
        }
    }
}
