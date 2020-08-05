#include "Functional.hxx"
#include "get_new_points.hxx"
#include "solve_LP.hxx"

bool is_feasible(const Functional &functional)
{
  size_t num_weights(functional.blocks.at(0).polys.size());

  const size_t num_blocks(functional.blocks.size());
  std::vector<El::BigFloat> weights(num_weights, 1);

  std::vector<std::set<El::BigFloat>> points(num_blocks);
  std::vector<std::vector<El::BigFloat>> new_points(num_blocks);

  const El::BigFloat scalar_gap(1.44);
  const int64_t max_twist(50), spacetime_dim(3);
  points.at(0).emplace(0);
  points.at(0).emplace(scalar_gap);
  new_points.at(0).emplace_back(max_twist);
  for(size_t block(1); block < num_blocks; ++block)
    {
      points.at(block).emplace(2 * block + spacetime_dim - 2);
      new_points.at(block).emplace_back(2 * block + max_twist);
    }

  const El::BigFloat min_weight_size(
    Sqrt(Sqrt(El::limits::Epsilon<El::BigFloat>())));
  bool has_new_points(true), is_too_small(false);
  while(has_new_points && !is_too_small)
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
          // std::cout << "points: " << block << " " << points.at(block) << "\n";
        }

      std::cout << "num_constraints: " << num_constraints << "\n";
      const size_t num_rows(num_constraints),
        num_columns(2 * weights.size() + num_constraints + 1);

      El::Matrix<El::BigFloat> A(num_rows, num_columns);
      El::Zero(A);

      size_t row(0);
      for(size_t block(0); block != num_blocks; ++block)
        {
          // One constraint per point
          for(auto &x : points.at(block))
            {
              // delta term
              A(row, 2 * weights.size()) = 1;
              // slack term
              A(row, 2 * weights.size() + row + 1) = -1;

              const El::BigFloat prefactor(
                functional.prefactor(x)
                * functional.blocks.at(block).pole_prefactor(x));
              for(size_t poly_index(0);
                  poly_index != functional.blocks.at(block).polys.size();
                  ++poly_index)
                {
                  A(row, 2 * poly_index)
                    = prefactor
                      * functional.blocks.at(block).polys.at(poly_index)(x);
                  A(row, 2 * poly_index + 1) = -A(row, 2 * poly_index);
                }
              ++row;
            }
        }

      El::Matrix<El::BigFloat> b(num_rows, 1), c(num_columns, 1);
      El::Zero(b);
      El::Zero(c);
      c(2 * weights.size(), 0) = 1;

      solve_LP(A, b, c, weights);

      // std::cout << "weight: " << weights << "\n";
      for(size_t block(0); block != num_blocks; ++block)
        {
          const El::BigFloat min_x(block == 0 ? scalar_gap
                                              : *(points.at(block).begin()));
          // 0.01 should be a small enough relative error so that we are
          // in the regime of convergence.  Then the error estimates will
          // work
          Mesh mesh(
            min_x, *(points.at(block).rbegin()),
            // *(points.at(block).begin()), *(points.at(block).rbegin()),
            [&](const El::BigFloat &x) {
              return functional.prefactor(x)
                     * functional.blocks.at(block).eval_weighted(x, weights);
            },
            0.01);
          new_points.at(block) = get_new_points(mesh);
          // std::cout << "new: " << block << " " << new_points.at(block) <<
          // "\n";
          has_new_points = has_new_points || !new_points.at(block).empty();
        }

      is_too_small = true;
      for(auto &weight : weights)
        {
          is_too_small = is_too_small && Abs(weight) < min_weight_size;
        }
    }
  // std::cout.precision(precision / 3.3);
  std::cout << "weights: " << weights << "\n";
  return !is_too_small;
}
