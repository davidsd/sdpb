#include "Mesh.hxx"
#include "Functional.hxx"

#include "../ostream_vector.hxx"

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
  const int64_t precision(1024);
  El::gmp::SetPrecision(precision);

  Functional functional("test/toy_polys", "", false);
  size_t num_weights(functional.blocks.at(0).polys.size());

  const size_t num_blocks(functional.blocks.size());
  const El::BigFloat min_x(0.0), max_x(100.0);
  std::vector<El::BigFloat> scale(num_blocks);
  std::vector<El::BigFloat> weights(num_weights, 1);
  for(size_t block(0); block != num_blocks; ++block)
    {
      scale[block]
        = Abs(functional.blocks[block].eval_weighted(max_x, weights));
    }

  std::vector<std::vector<El::BigFloat>> points(num_blocks, {min_x}),
    new_points(num_blocks, {max_x});

  bool has_new_points(true);
  while(has_new_points)
    {
      has_new_points = false;
      size_t num_constraints(0);
      for(size_t block(0); block != num_blocks; ++block)
        {
          // 0.01 is completely arbitrary.  We want it big enough to not
          // rule out points that might provide a limit, but not so big to
          // not rule out any points.

          // For the toy example, this eliminates almost all of the
          // points.  It feels dangerous.
          El::BigFloat tolerance(
            scale.at(block) * Sqrt(Sqrt(El::limits::Epsilon<El::BigFloat>())));
          {
            std::vector<El::BigFloat> temp_points;
            El::BigFloat f0, f1;
            for(auto &point : points.at(block))
              {
                // if(functional.blocks.at(block).eval_weighted(point, weights)
                //    < tolerance)
                  {
                    temp_points.emplace_back(point);
                  }
              }
            std::swap(points.at(block), temp_points);
          }
          for(auto &point : new_points.at(block))
            {
              points.at(block).emplace_back(point);
            }
          num_constraints += points.at(block).size();
        }
      
      const size_t num_rows(num_constraints + 1),
        // const size_t num_rows(num_constraints),
        num_columns(2 * weights.size() + num_constraints + 1);

      El::Matrix<El::BigFloat> A(num_rows, num_columns);
      El::Zero(A);

      size_t row(0);
      for(size_t block(0); block != num_blocks; ++block)
        {
          // One constraint per point
          for(size_t point(0); point != points.at(block).size(); ++point)
            {
              // delta term
              A(row, 0) = 1;
              // slack term
              A(row, 2 * weights.size() + row + 1) = -1;

              const El::BigFloat &x(points.at(block)[point]);
              for(size_t poly_index(0);
                  poly_index != functional.blocks.at(block).polys.size();
                  ++poly_index)
                {
                  A(row, 2 * poly_index + 1)
                    = functional.blocks.at(block).polys.at(poly_index)(x);
                  A(row, 2 * poly_index + 2) = -A(row, 2 * poly_index + 1);
                }
              ++row;
            }
        }
      A(num_rows - 1, 1) = 1;
      A(num_rows - 1, 2) = 1;

      El::Matrix<El::BigFloat> b(num_rows, 1), c(num_columns, 1);
      El::Zero(b);
      b(num_rows - 1, 0) = 1;
      El::Zero(c);
      c(0, 0) = 1;

      // std::cout << "A: "
      //           << A.Height() << " "
      //           << A.Width() << " "
      //           << b.Height() << " "
      //           << b.Width() << " "
      //           << c.Height() << " "
      //           << c.Width() << " "
      //           << "\n";

      // El::Print(A, "A");
      // El::Print(b, "\nb");
      // El::Print(c, "\nc");
      // std::cout << "\n";

      solve_LP(A, b, c, weights);

      std::cout << "weight: " << weights << "\n";
      for(size_t block(0); block != num_blocks; ++block)
        {
          // std::cout.precision(precision / 3.3);
          // std::cout << "solve: " << points.at(block).size() << " " << weights
          //           << "\n";

          // 0.01 should be a small enough relative error so that we are
          // in the regime of convergence.  Then the error estimates will
          // work
          Mesh mesh(min_x, max_x,
                    [&](const El::BigFloat &x) {
                      return functional.blocks.at(block).eval_weighted(
                        x, weights);
                    },
                    0.01);
          new_points.at(block) = get_new_points(mesh);
          std::cout << "new: " << new_points.at(block) << "\n";
          scale.at(block) = max_value(mesh);
          has_new_points = has_new_points || !new_points.at(block).empty();
        }
    }
  // std::cout.precision(precision / 3.3);
  std::cout << "weights: " << weights << "\n";
}
