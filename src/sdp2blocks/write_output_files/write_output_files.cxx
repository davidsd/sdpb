#include "../Positive_Matrix_With_Prefactor.hxx"

#include <boost/multiprecision/gmp.hpp>
#include <boost/filesystem.hpp>
#include <algorithm>

using Boost_Float = boost::multiprecision::mpfr_float;
std::vector<Boost_Float> sample_points(const size_t &num_points);

void write_output_files(
  const boost::filesystem::path &outdir,
  const std::vector<El::BigFloat> &objectives,
  const std::vector<El::BigFloat> &normalization,
  const std::vector<Positive_Matrix_With_Prefactor> &matrices)
{
  auto max_normalization(
    std::max_element(normalization.begin(), normalization.end()));
  size_t max_index(std::distance(normalization.begin(), max_normalization));

  El::BigFloat objective_const(objectives.at(max_index)
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

  for(auto &matrix : matrices)
    {
      auto max_polynomial(
        std::max_element(matrix.polynomials.begin(), matrix.polynomials.end(),
                         [](const Polynomial &a, const Polynomial &b) {
                           return a.degree() < b.degree();
                         }));
      if(max_polynomial == matrix.polynomials.end())
        {
          throw std::runtime_error("No polynomials in the Positive Matrix");
        }

      const size_t max_degree(max_polynomial->degree());
      size_t rows(matrix.polynomials.size()), cols(max_degree + 1);
      std::vector<Boost_Float> points(sample_points(max_degree + 1)),
        sample_scalings;

      std::cout << "Points: " << (max_degree + 1) << " " << points[0] << "\n";

      sample_scalings.reserve(points.size());
      for(auto &point : points)
        {
          Boost_Float numerator(matrix.damped_rational.constant
                                * pow(matrix.damped_rational.base, point));
          Boost_Float denominator(1);
          std::cout << "Const: " << matrix.damped_rational.constant << " "
                    << matrix.damped_rational.base << " " << point << "\n";
          for(auto &pole : matrix.damped_rational.poles)
            {
              denominator *= (point - pole);
            }
          sample_scalings.push_back(numerator / denominator);
          std::cout << "Scaling: " << sample_scalings.back() << "\n";
        }
    }

  // std::cout << "Objective: " << objective_const << "\n";
  // for(auto &o : dual_objective_b)
  //   std::cout << " " << o << "\n";
}
