#include "Mesh.hxx"
#include "Functional.hxx"

#include "../ostream_sequence.hxx"

// We convert the optimization problem into a regular linear
// programming problem.
//
// 1) Each polynomial in Block::polys adds two variables.  The
// weights of those polynomials are unbounded, while linear
// programming requires a strictly positive variable.  So we
// substitute W_n = w_n+ - w_n+, where both w_n+ and w_n- are
// strictly positive.  All blocks should have the same number of
// polynomials, so we only need to look at the first one.
//
// 2) Each constraint adds one 'slack' variable s_n.  There is one
// constraint per coordinate 'x', with each block having multiple,
// separate coordinates.
//
// 3) One more global variable delta, which gives the linear
// program something to minimize.
//
// This turns the problem
//
//   A_0 . W > 0
//   A_1 . W > 0
//   ...
//
// into
//
//   min delta
//
// where
//
//   A_0 . (w_+ - w_-) + delta - s_0 = 0
//   A_1 . (w_+ - w_-) + delta - s_1 = 0
//   ...
//   w_n+, w_n-, s_n, delta >= 0
//
// There is a constraint for every point that is sampled.  At the
// beginning, we sample the min and max for each block, so there are
// 2*num_blocks constraints.  In general
//
// => num_rows == num_constraints
// => num_columns == 2*num_weights + num_constraints + 1

std::vector<El::BigFloat> get_new_points(const Mesh &mesh);
El::BigFloat max_value(const Mesh &mesh);

void solve_LP(const El::Matrix<El::BigFloat> &A,
              const El::Matrix<El::BigFloat> &b,
              const El::Matrix<El::BigFloat> &c,
              std::vector<El::BigFloat> &weights);

int main(int argc, char **argv)
{
  El::Environment env(argc, argv);
  const int64_t precision(256);
  El::gmp::SetPrecision(precision);

  Functional functional("test/single_corr_polys", "test/single_corr_poles");
  // Functional functional("test/toy_polys", "");
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

  bool has_new_points(true);
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

      const size_t num_rows(num_constraints + 1),
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
              A(row, 0) = 1;
              // slack term
              A(row, 2 * weights.size() + row + 1) = -1;

              const El::BigFloat prefactor(
                functional.prefactor(x)
                * functional.blocks.at(block).pole_prefactor(x));
              for(size_t poly_index(0);
                  poly_index != functional.blocks.at(block).polys.size();
                  ++poly_index)
                {
                  A(row, 2 * poly_index + 1)
                    = prefactor
                      * functional.blocks.at(block).polys.at(poly_index)(x);
                  A(row, 2 * poly_index + 2) = -A(row, 2 * poly_index + 1);
                }
              ++row;
            }
        }

      El::Matrix<El::BigFloat> b(num_rows, 1), c(num_columns, 1);
      El::Zero(b);
      El::Zero(c);
      c(0, 0) = 1;

      // Set the first coefficient to 1
      A(num_rows - 1, 1) = 1;
      A(num_rows - 1, 2) = -1;
      b(num_rows - 1) = 1;

      solve_LP(A, b, c, weights);

      std::cout << "weight: " << weights << "\n";
      for(size_t block(0); block != num_blocks; ++block)
        {
          const El::BigFloat min_x(block==0 ? scalar_gap : *(points.at(block).begin()));
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
    }
  // std::cout.precision(precision / 3.3);
  std::cout << "weights: " << weights << "\n";
}
