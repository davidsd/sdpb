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
  El::gmp::SetPrecision(256);

  El::Matrix<El::BigFloat> b(1, 1);
  El::Matrix<El::BigFloat> A(1, 3), c(3, 1);

  c(0) = 1;
  c(1) = -1;

  std::vector<El::BigFloat> points({0.0});
  eval(points[0], b(0), A(0, 1));
  A(0, 0) = -A(0, 1);
  A(0, 2) = 1;

  std::vector<El::BigFloat> new_points({100.0});

  while(!new_points.empty())
    {
      const size_t old_size(points.size()), delta(new_points.size()),
        new_size(old_size + delta);
      points.reserve(new_size);
      {
        El::Matrix<El::BigFloat> b_new(new_size, 1);
        El::Matrix<El::BigFloat> A_new(new_size, new_size + 2),
          c_new(new_size + 2, 1);

        for(size_t index_0(0); index_0 < old_size; ++index_0)
          {
            b_new(index_0) = b(index_0);
          }
        for(size_t index_0(0); index_0 < old_size; ++index_0)
          for(size_t index_1(0); index_1 < old_size + 2; ++index_1)
            {
              A_new(index_0, index_1) = A(index_0, index_1);
            }
        for(size_t index_1(0); index_1 < old_size + 2; ++index_1)
          {
            c_new(index_1) = c(index_1);
          }
        A = A_new;
        b = b_new;
        c = c_new;
      }
      for(size_t d(0); d != delta; ++d)
        {
          c(old_size + d + 2) = 0;
          for(size_t index(0); index != old_size; ++index)
            {
              A(index, old_size + d + 2) = 0;
            }

          eval(new_points[d], b(old_size + d), A(old_size + d, 1));
          A(old_size + d, 0) = -A(old_size + d, 1);

          for(size_t index(2); index != old_size + d + 2; ++index)
            {
              A(old_size + d, index) = 0;
            }
          A(old_size + d, old_size + d + 2) = 1;
          points.emplace_back(new_points[d]);
        }

      // El::Print(A, "A");
      // std::cout << "\n";
      // El::Print(b, "b");
      // std::cout << "\n";
      // El::Print(c, "c");
      // std::cout << "\n";

      // std::cout << "solve: " << solve_LP(A, b, c) << "\n";

      El::BigFloat optimal(solve_LP(A, b, c));
      std::cout << "solve: " << optimal << "\n";
      Mesh mesh(points[0],points[1],[=](const El::BigFloat &x)
                          {
                            El::BigFloat f0,f1;
                            eval(x,f0,f1);
                            return f0 + optimal * f1;
                          },
        0.01);
      // std::cout << "mesh: "
      //           << mesh << "\n";

      new_points = get_new_points(mesh);

      std::cout << "new: " << new_points << "\n";
    }
}
