#include "Mesh.hxx"

#include "../ostream_vector.hxx"

void eval(const El::BigFloat &x, El::BigFloat &f0, El::BigFloat &f1)
{
  const El::BigFloat pow2(x * x), pow4(pow2 * pow2);
  f0 = 1 + pow4;
  f1 = pow2 + pow4 / 12;
}

std::vector<El::BigFloat> get_new_points(const Mesh &mesh);

El::BigFloat
solve_LP(const El::Matrix<El::BigFloat> &A, const El::Matrix<El::BigFloat> &b,
         const El::Matrix<El::BigFloat> &c);

int main(int argc, char **argv)
{
  El::Environment env(argc, argv);
  El::gmp::SetPrecision(1024);

  std::vector<El::BigFloat> points({0.0}), new_points({100.0});

  while(!new_points.empty())
    {
      for(auto &point : new_points)
        {
          points.emplace_back(point);
        }
      El::Matrix<El::BigFloat> b(points.size(), 1),
        A(points.size(), points.size() + 2), c(points.size() + 2, 1);
      c(0) = 1;
      c(1) = -1;

      for(size_t point(0); point != points.size(); ++point)
        {
          c(point + 2) = 0;
          for(size_t index(0); index != points.size(); ++index)
            {
              A(index, point + 2) = 0;
            }
          A(point, point + 2) = 1;

          eval(points[point], b(point), A(point, 1));
          A(point, 0) = -A(point, 1);
        }
      El::BigFloat optimal(solve_LP(A, b, c));
      std::cout.precision(1024 / 3.3);
      std::cout << "solve: " << optimal << "\n";
      Mesh mesh(points[0], points[1],
                [=](const El::BigFloat &x) {
                  El::BigFloat f0, f1;
                  eval(x, f0, f1);
                  return f0 + optimal * f1;
                },
                0.01);
      new_points = get_new_points(mesh);
      std::cout << "new: " << new_points << "\n";
    }
}
