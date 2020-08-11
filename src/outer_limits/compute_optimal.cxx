#include "get_new_points.hxx"
#include "solve_LP.hxx"

#include "eval_weighted.hxx"
#include "poles_prefactor.hxx"
#include "power_prefactor.hxx"

#include "../sdp_solve.hxx"
#include "../ostream_set.hxx"

std::vector<El::BigFloat>
compute_optimal(const std::vector<Positive_Matrix_With_Prefactor> &matrices,
                const std::vector<El::BigFloat> &objectives,
                const std::vector<El::BigFloat> &normalization)
{
  size_t num_weights(normalization.size());

  const size_t num_blocks(matrices.size());
  std::vector<El::BigFloat> weights(num_weights, 1);

  std::vector<std::set<El::BigFloat>> points(num_blocks);
  std::vector<std::vector<El::BigFloat>> new_points(num_blocks);

  std::cout << "blocks: " << num_blocks << "\n";
  const El::BigFloat min_x(0), max_x(10);
  for(size_t block(0); block < num_blocks; ++block)
    {
      points.at(block).emplace(min_x);
      new_points.at(block).emplace_back(max_x);
    }

  bool has_new_points(true);

  const size_t procs_per_node(1), proc_granularity(1);
  const Verbosity verbosity(Verbosity::regular);
  while(has_new_points)
    {
      has_new_points = false;
      size_t num_constraints(0);
      for(size_t block(0); block != num_blocks; ++block)
        {
          for(auto &point : new_points.at(block))
            {
              points.at(block).emplace(point);
            }
          num_constraints += points.at(block).size();
          std::cout << "points: " << block << " " << points.at(block) << "\n";
        }

      std::cout << "num_constraints: " << num_constraints << "\n";

      Block_Info block_info(num_constraints, procs_per_node, proc_granularity,
                            verbosity);
      El::Grid grid(block_info.mpi_comm.value);
      std::vector<El::BigFloat> prefactor;
      for(size_t block(0); block != num_blocks; ++block)
        {
          for(auto &point : points.at(block))
            {
              prefactor.push_back(
                power_prefactor(matrices[block].damped_rational.base, point)
                * poles_prefactor(matrices[block].damped_rational.poles, point));
            }
        }
      SDP sdp(objectives, normalization, prefactor, block_info, grid);

      const size_t num_rows(num_constraints + 1),
        num_columns(2 * weights.size() + num_constraints);

      El::Matrix<El::BigFloat> A(num_rows, num_columns);
      El::Zero(A);

      size_t row(0);
      for(size_t block(0); block != num_blocks; ++block)
        {
          // One constraint per point
          for(auto &x : points.at(block))
            {
              const El::BigFloat prefactor(
                power_prefactor(matrices[block].damped_rational.base, x)
                * poles_prefactor(matrices[block].damped_rational.poles, x));

              for(size_t matrix_row(0);
                  matrix_row < matrices[block].polynomials.size();
                  ++matrix_row)
                for(size_t matrix_column(0);
                    matrix_column
                    < matrices[block].polynomials[matrix_row].size();
                    ++matrix_column)
                  {
                    // slack term
                    A(row, 2 * weights.size() + row) = -1;

                    for(size_t poly_index(0);
                        poly_index
                        != matrices[block]
                             .polynomials[matrix_row][matrix_column]
                             .size();
                        ++poly_index)
                      {
                        A(row, 2 * poly_index)
                          = prefactor
                            * matrices[block]
                                .polynomials[matrix_row][matrix_column]
                                            [poly_index](x);
                        A(row, 2 * poly_index + 1) = -A(row, 2 * poly_index);
                      }
                    ++row;
                  }
            }
        }

      for(size_t index(0); index != normalization.size(); ++index)
        {
          A(num_rows - 1, 2 * index) = normalization[index];
          A(num_rows - 1, 2 * index + 1) = -normalization[index];
        }

      El::Matrix<El::BigFloat> b(num_rows, 1);
      El::Zero(b);
      b(num_rows - 1, 0) = 1;

      El::Matrix<El::BigFloat> c(num_columns, 1);
      El::Zero(c);
      for(size_t index(0); index != objectives.size(); ++index)
        {
          c(2 * index, 0) = -objectives[index];
          c(2 * index + 1, 0) = objectives[index];
        }

      solve_LP(A, b, c, weights);

      // std::cout << "weight: " << weights << "\n";
      for(size_t block(0); block != num_blocks; ++block)
        {
          // 0.01 should be a small enough relative error so that we are
          // in the regime of convergence.  Then the error estimates will
          // work
          Mesh mesh(*(points.at(block).begin()), *(points.at(block).rbegin()),
                    [&](const El::BigFloat &x) {
                      return eval_weighted(matrices[block], x, weights);
                    },
                    0.01);
          new_points.at(block) = get_new_points(mesh);
          // std::cout << "new: " << block << " " << new_points.at(block) <<
          // "\n";
          has_new_points = has_new_points || !new_points.at(block).empty();
        }
    }
  // std::cout.precision(precision / 3.3);
  std::cout << "weights: " << weights << "\n";
  return weights;
}
