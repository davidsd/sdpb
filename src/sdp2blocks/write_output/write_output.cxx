#include "../Positive_Matrix_With_Prefactor.hxx"
#include "../Boost_Float.hxx"
#include "../../Timers.hxx"

#include <boost/filesystem.hpp>
#include <algorithm>

El::Matrix<El::BigFloat>
bilinear_basis(const Damped_Rational &damped_rational,
               const size_t &half_max_degree, const std::string &timer_prefix,
               Timers &timers);

std::vector<Boost_Float> sample_points(const size_t &num_points);

void write_output(const boost::filesystem::path &outdir,
                  const std::vector<El::BigFloat> &objectives,
                  const std::vector<El::BigFloat> &normalization,
                  const std::vector<Positive_Matrix_With_Prefactor> &matrices,
                  Timers &timers)
{
  auto &objectives_timer(timers.add_and_start("write_output.objectives"));
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

  objectives_timer.stop();

  auto &matrices_timer(timers.add_and_start("write_output.matrices"));
  for(size_t index = 0; index != matrices.size(); ++index)
    {
      auto &scalings_timer(timers.add_and_start(
        "write_output.matrices.scalings_" + std::to_string(index)));
      const size_t max_degree([&]() {
        int64_t result(0);
        for(auto &polynomial : matrices[index].polynomials)
          {
            result = std::max(result, polynomial.degree());
          }
        return result;
      }());
      size_t rows(matrices[index].polynomials.size()), cols(max_degree + 1);
      std::vector<Boost_Float> points(sample_points(max_degree + 1)),
        sample_scalings;

      sample_scalings.reserve(points.size());
      for(auto &point : points)
        {
          Boost_Float numerator(
            matrices[index].damped_rational.constant
            * pow(matrices[index].damped_rational.base, point));
          Boost_Float denominator(1);
          for(auto &pole : matrices[index].damped_rational.poles)
            {
              denominator *= (point - pole);
            }
          sample_scalings.push_back(numerator / denominator);
        }
      scalings_timer.stop();

      const std::string bilinear_basis_timer_name(
        "write_output.matrices.bilinear_basis_" + std::to_string(index));
      auto &bilinear_basis_timer(
        timers.add_and_start(bilinear_basis_timer_name));
      bilinear_basis(matrices[index].damped_rational, max_degree / 2,
                     bilinear_basis_timer_name, timers);
      bilinear_basis_timer.stop();
    }
  matrices_timer.stop();

  // std::cout << "Objective: " << objective_const << "\n";
  // for(auto &o : dual_objective_b)
  //   std::cout << " " << o << "\n";
}
