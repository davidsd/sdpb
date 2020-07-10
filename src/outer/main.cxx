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

  const El::BigFloat min_x(0.0), max_x(100.0);
  std::vector<El::BigFloat> points({min_x}), new_points({max_x});

  El::BigFloat f0, f1, scale;
  eval(min_x, f0, f1);
  scale = Max(Abs(f0), Abs(f1));
  eval(max_x, f0, f1);
  scale = Max(scale, Max(Abs(f0), Abs(f1)));

  // {
  //   Functional functional("test/toy_polys","",false);

  //   std::vector<std::vector<El::BigFloat>> optimals(1);
  //   optimals[0].resize(2);
  //   optimals[0][0]=1;
  //   optimals[0][1]=1;

  //   std::vector<El::BigFloat> coords(1);
  //   coords[0]=min_x;
  //   std::vector<El::BigFloat> result(functional.eval(coords,optimals));
  //   std::cout << "functional: " << result[0] << "\n";
  //   coords[0]=max_x;
  //   std::cout << "functional: " << functional.eval(coords,optimals)[0] << "\n";
  //   exit(0);
  // }
  
  El::BigFloat optimal(0);
  while(!new_points.empty())
    {
      // 0.01 is completely arbitrary.  We want it big enough to not
      // rule out points that might provide a limit, but not so big to
      // not rule out any points.

      // For the toy example, this eliminates almost all of the
      // points.  It feels dangerous.
      El::BigFloat tolerance(
        scale * Sqrt(Sqrt(El::limits::Epsilon<El::BigFloat>())));
      {
        std::vector<El::BigFloat> temp_points;
        El::BigFloat f0, f1;
        for(auto &point : points)
          {
            eval(point, f0, f1);
            if(f0 + optimal * f1 < tolerance)
              {
                temp_points.emplace_back(point);
              }
          }
        std::swap(points, temp_points);
      }
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
      optimal = solve_LP(A, b, c);
      std::cout.precision(precision / 3.3);
      std::cout << "solve: " << points.size() << " " << optimal << "\n";
      // 0.01 should be a small enough relative error so that we are
      // in the regime of convergence.  Then the error estimates will
      // work
      Mesh mesh(min_x, max_x,
                [=](const El::BigFloat &x) {
                  El::BigFloat f0, f1;
                  eval(x, f0, f1);
                  return f0 + optimal * f1;
                },
                0.01);
      new_points = get_new_points(mesh);
      // std::cout << "new: " << new_points << "\n";
      scale = max_value(mesh);
    }
  std::cout.precision(precision / 3.3);
  std::cout << "optimal: " << optimal << "\n";
}
