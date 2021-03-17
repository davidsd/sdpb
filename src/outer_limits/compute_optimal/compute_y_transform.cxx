#include "Mesh.hxx"

#include "../Function.hxx"
#include "../../sdp_solve.hxx"
#include "../../ostream_set.hxx"
#include "../../ostream_vector.hxx"
#include "../../set_stream_precision.hxx"

void setup_constraints(
  const size_t &max_index, const size_t &num_blocks,
  const El::BigFloat &infinity,
  const std::vector<std::vector<std::vector<std::vector<Function>>>>
    &function_blocks,
  const std::vector<El::BigFloat> &normalization,
  const std::vector<std::set<El::BigFloat>> &points,
  std::vector<std::vector<El::BigFloat>> &primal_objective_c,
  std::vector<El::Matrix<El::BigFloat>> &free_var_matrix);

void compute_y_transform(
  const std::vector<std::vector<std::vector<std::vector<Function>>>>
    &function_blocks,
  const std::vector<std::set<El::BigFloat>> &points,
  const std::vector<El::BigFloat> &objectives,
  const std::vector<El::BigFloat> &normalization,
  const SDP_Solver_Parameters &parameters, const size_t &max_index,
  const El::Grid &global_grid,
  El::DistMatrix<El::BigFloat, El::STAR, El::STAR> &yp_to_y_star,
  El::DistMatrix<El::BigFloat, El::STAR, El::STAR> &dual_objective_b_star,
  El::BigFloat &primal_c_scale)
{
  const size_t num_blocks(points.size());

  // GMP does not have a special infinity value, so we use max double.
  const El::BigFloat infinity(std::numeric_limits<double>::max());

  size_t num_constraints(0);
  std::vector<size_t> matrix_dimensions;
  for(size_t block(0); block != num_blocks; ++block)
    {
      num_constraints += points.at(block).size();
      matrix_dimensions.insert(matrix_dimensions.end(),
                               points.at(block).size(),
                               function_blocks[block].size());
    }

  Block_Info block_info(matrix_dimensions, parameters.procs_per_node,
                        parameters.proc_granularity, parameters.verbosity);

  std::vector<std::vector<El::BigFloat>> primal_objective_c;
  primal_objective_c.reserve(num_constraints);
  std::vector<El::Matrix<El::BigFloat>> free_var_matrix;
  free_var_matrix.reserve(num_constraints);

  setup_constraints(max_index, num_blocks, infinity, function_blocks,
                    normalization, points, primal_objective_c,
                    free_var_matrix);

  const El::BigFloat objective_const(objectives.at(max_index)
                                     / normalization.at(max_index));
  std::vector<El::BigFloat> dual_objective_b;
  dual_objective_b.reserve(normalization.size() - 1);
  for(size_t index = 0; index < normalization.size(); ++index)
    {
      if(index != max_index)
        {
          dual_objective_b.push_back(
            objectives.at(index) - normalization.at(index) * objective_const);
        }
    }

  auto &block_indices(block_info.block_indices);
  El::BigFloat max_c(0.0);
  for(size_t block(0); block != block_indices.size(); ++block)
    {
      for(auto &element :
          primal_objective_c.at(block_info.block_indices.at(block)))
        {
          max_c = std::max(max_c, El::Abs(element));
        }
    }
  max_c = El::mpi::AllReduce(max_c, El::mpi::MAX, El::mpi::COMM_WORLD);
  primal_c_scale
    = (max_c == El::BigFloat(0.0)) ? El::BigFloat(1.0) : 1 / max_c;

  // Setup B
  std::vector<size_t> block_offsets(primal_objective_c.size() + 1, 0);
  for(size_t p(0); p < primal_objective_c.size(); ++p)
    {
      block_offsets[p + 1] = block_offsets[p] + primal_objective_c[p].size();
    }

  const int64_t B_Height(block_offsets.back()),
    B_Width(dual_objective_b.size());
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
              * free_var_matrix.at(block_index)(block_row, global_column));
        }
    }

  // Compute SVD of B
  //   B = U s V^T
  // Define
  //   b''(j) = Sum(b(m) * V^T(j,m)/s(j), m)
  //   B'' = U
  // This gives
  //   minimize b''(j) y''(j) + c_0
  // with constraints
  //   B''(i,j) y''(j) <= c(i)
  // This implies
  //   y(m) = Sum((V^T(l,m)/s(l)) * y''(l), l)
  // and so to convert back to y
  //   y(m) = Sum(yp_to_y(m,l) * y''(l), l)
  //   yp_to_y(m,l) = V^T(l,m)/s(l)

  El::DistMatrix<El::BigFloat> U(global_grid);
  El::DistMatrix<El::BigFloat> temp(global_grid), V(global_grid),
    dual_objective_b_global(dual_objective_b.size(), 1, global_grid);
  // SVD returns U, s, and V
  El::SVD(B, U, temp, V);

  El::DistMatrix<El::BigFloat> yp_to_y(V);
  El::DiagonalSolve(El::LeftOrRight::RIGHT, El::Orientation::NORMAL, temp,
                    yp_to_y);
  yp_to_y_star = yp_to_y;

  for(int64_t row(0); row < temp.LocalHeight(); ++row)
    {
      const int64_t global_row(temp.GlobalRow(row));
      for(int64_t column(0); column < temp.LocalWidth(); ++column)
        {
          temp.SetLocal(row, column, dual_objective_b.at(global_row));
        }
    }

  El::Zero(dual_objective_b_global);
  El::Gemv(El::Orientation::TRANSPOSE, El::BigFloat(1.0), yp_to_y_star, temp,
           El::BigFloat(0.0), dual_objective_b_global);

  const El::BigFloat max_b(El::MaxAbs(dual_objective_b_global));
  const El::BigFloat b_scale((max_b == El::BigFloat(0.0)) ? El::BigFloat(1.0)
                                                          : 1 / max_b);

  dual_objective_b_global *= b_scale;
  dual_objective_b_star = dual_objective_b_global;
}
