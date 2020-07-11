#include "Mesh.hxx"
#include "Functional.hxx"

#include "../ostream_vector.hxx"

void eval(const El::BigFloat &x, El::BigFloat &f0, El::BigFloat &f1)
{
  const El::BigFloat pow2(x * x), pow4(pow2 * pow2);
  f0 = 1 + pow4;
  f1 = pow2 + pow4 / 12;
}

std::vector<El::BigFloat> get_new_points(const Mesh &mesh);
El::BigFloat max_value(const Mesh &mesh);

El::BigFloat
solve_LP(const El::Matrix<El::BigFloat> &A, const El::Matrix<El::BigFloat> &b,
         const El::Matrix<El::BigFloat> &c);

int main(int argc, char **argv)
{
  El::Environment env(argc, argv);
  const int64_t precision(1024);
  El::gmp::SetPrecision(precision);

  Functional functional("test/toy_polys", "", false);

  const El::BigFloat min_x(0.0), max_x(100.0);
  std::vector<El::BigFloat> scale(functional.blocks.size());
  std::vector<std::vector<El::BigFloat>> weights(functional.blocks.size());
  for(size_t index(0); index != functional.blocks.size(); ++index)
    {
      weights[index].resize(functional.blocks[index].p.size(), 1);
      scale[index] = Abs(functional.blocks[index].eval(max_x, weights[index]));
    }

  std::vector<std::vector<El::BigFloat>> points(scale.size(), {min_x}),
    new_points(scale.size(), {max_x});

  bool has_new_points(true);
  while(has_new_points)
    {
      has_new_points = false;
      for(size_t index(0); index != scale.size(); ++index)
        {
          // 0.01 is completely arbitrary.  We want it big enough to not
          // rule out points that might provide a limit, but not so big to
          // not rule out any points.

          // For the toy example, this eliminates almost all of the
          // points.  It feels dangerous.
          El::BigFloat tolerance(
            scale.at(index) * Sqrt(Sqrt(El::limits::Epsilon<El::BigFloat>())));
          {
            std::vector<El::BigFloat> temp_points;
            El::BigFloat f0, f1;
            for(auto &point : points.at(index))
              {
                if(functional.blocks.at(index).eval(point, weights.at(index))
                   < tolerance)
                  {
                    temp_points.emplace_back(point);
                  }
              }
            std::swap(points.at(index), temp_points);
          }
          for(auto &point : new_points.at(index))
            {
              points.at(index).emplace_back(point);
            }
          El::Matrix<El::BigFloat> b(points.at(index).size(), 1),
            A(points.at(index).size(), points.at(index).size() + 2),
            c(points.at(index).size() + 2, 1);
          c(0) = 1;
          c(1) = -1;

          for(size_t point(0); point != points.at(index).size(); ++point)
            {
              c(point + 2) = 0;
              for(size_t row(0); row != points.at(index).size(); ++row)
                {
                  A(row, point + 2) = 0;
                }
              A(point, point + 2) = 1;

              eval(points.at(index)[point], b(point), A(point, 1));
              A(point, 0) = -A(point, 1);
            }
          weights.at(index)[1] = solve_LP(A, b, c);
          std::cout.precision(precision / 3.3);
          std::cout << "solve: " << points.at(index).size() << " "
                    << weights.at(index)[1] << "\n";
          // 0.01 should be a small enough relative error so that we are
          // in the regime of convergence.  Then the error estimates will
          // work
          Mesh mesh(min_x, max_x,
                    [&](const El::BigFloat &x) {
                      return functional.blocks.at(index).eval(
                        x, weights.at(index));
                    },
                    0.01);
          new_points.at(index) = get_new_points(mesh);
          // std::cout << "new: " << new_points.at(index) << "\n";
          scale.at(index) = max_value(mesh);
          has_new_points = has_new_points || !new_points.at(index).empty();
        }
    }
  std::cout.precision(precision / 3.3);
  std::cout << "weights: " << weights << "\n";
}
