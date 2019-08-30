#include "../Positive_Matrix_With_Prefactor.hxx"
#include "../Boost_Float.hxx"
#include "../../Timers.hxx"
#include "../../sdp_convert.hxx"

#include <boost/filesystem.hpp>
#include <algorithm>

std::vector<Polynomial> bilinear_basis(const Damped_Rational &damped_rational,
                                       const size_t &half_max_degree);

std::vector<Boost_Float> sample_points(const size_t &num_points);

void write_output(const boost::filesystem::path &output_dir,
                  const std::vector<El::BigFloat> &objectives,
                  const std::vector<El::BigFloat> &normalization,
                  const std::vector<Positive_Matrix_With_Prefactor> &matrices,
                  Timers &timers)
{
  auto &objectives_timer(timers.add_and_start("write_output.objectives"));

  auto max_normalization(normalization.begin());
  for(auto n(normalization.begin()); n!=normalization.end(); ++n)
    {
      if(Abs(*n)>Abs(*max_normalization))
        {
          max_normalization=n;
        }
    }
  size_t max_index(std::distance(normalization.begin(), max_normalization));

  El::BigFloat objective_const(objectives.at(max_index)
                               / normalization.at(max_index));
  std::vector<El::BigFloat> dual_objectives_b;
  dual_objectives_b.reserve(normalization.size() - 1);
  for(size_t index = 0; index < normalization.size(); ++index)
    {
      if(index != max_index)
        {
          dual_objectives_b.push_back(
            objectives.at(index) - normalization.at(index) * objective_const);
        }
    }

  objectives_timer.stop();

  auto &matrices_timer(timers.add_and_start("write_output.matrices"));
  std::vector<Dual_Constraint_Group> dual_constraint_groups;
  int rank(El::mpi::Rank(El::mpi::COMM_WORLD)),
    num_procs(El::mpi::Size(El::mpi::COMM_WORLD));
  std::vector<size_t> indices;
  for(size_t index = rank; index < matrices.size(); index += num_procs)
    {
      indices.push_back(index);
    }
  for(auto &index : indices)
    {
      auto &scalings_timer(timers.add_and_start(
        "write_output.matrices.scalings_" + std::to_string(index)));
      const size_t max_degree([&]() {
        int64_t result(0);
        for(auto &pvv : matrices[index].polynomials)
          for(auto &pv : pvv)
            for(auto &polynomial : pv)
              {
                result = std::max(result, polynomial.degree());
              }
        return result;
      }());
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

      Polynomial_Vector_Matrix pvm;
      pvm.rows = matrices[index].polynomials.size();
      pvm.cols = matrices[index].polynomials.front().size();

      auto &bilinear_basis_timer(timers.add_and_start(
        "write_output.matrices.bilinear_basis_" + std::to_string(index)));

      pvm.bilinear_basis
        = bilinear_basis(matrices[index].damped_rational, max_degree / 2);

      bilinear_basis_timer.stop();

      pvm.sample_points.reserve(points.size());
      for(auto &point : points)
        {
          pvm.sample_points.emplace_back(to_string(point));
        }
      pvm.sample_scalings.reserve(sample_scalings.size());
      for(auto &scaling : sample_scalings)
        {
          pvm.sample_scalings.emplace_back(to_string(scaling));
        }

      auto &pvm_timer(timers.add_and_start("write_output.matrices.pvm_"
                                           + std::to_string(index)));
      pvm.elements.reserve(pvm.rows * pvm.cols);
      for(auto &pvv : matrices[index].polynomials)
        for(auto &pv : pvv)
          {
            pvm.elements.emplace_back();
            auto &pvm_polynomials(pvm.elements.back());
            pvm_polynomials.reserve(pv.size());
            pvm_polynomials.push_back(pv.at(max_index)
                                      / normalization.at(max_index));
            auto &pvm_constant(pvm_polynomials.back());

            for(size_t index = 0; index < normalization.size(); ++index)
              {
                if(index != max_index)
                  {
                    pvm_polynomials.emplace_back(0, 0);
                    auto &pvm_poly(pvm_polynomials.back());
                    pvm_poly.coefficients.reserve(pv.at(index).degree() + 1);
                    size_t coefficient(0);
                    for(; coefficient < pv.at(index).coefficients.size()
                          && coefficient < pvm_constant.coefficients.size();
                        ++coefficient)
                      {
                        pvm_poly.coefficients.push_back(
                          pv.at(index).coefficients[coefficient]
                          - normalization.at(index)
                              * pvm_constant.coefficients[coefficient]);
                      }
                    for(; coefficient < pv.at(index).coefficients.size();
                        ++coefficient)
                      {
                        pvm_poly.coefficients.push_back(
                          pv.at(index).coefficients[coefficient]);
                      }
                    for(; coefficient < pvm_constant.coefficients.size();
                        ++coefficient)
                      {
                        pvm_poly.coefficients.push_back(
                          -normalization.at(index)
                          * pvm_polynomials.at(0).coefficients[coefficient]);
                      }
                  }
              }
          }
      pvm_timer.stop();
      auto &dual_constraint_timer(timers.add_and_start(
        "write_output.matrices.dual_constraint_" + std::to_string(index)));
      dual_constraint_groups.emplace_back(pvm);
      dual_constraint_timer.stop();
    }
  matrices_timer.stop();

  auto &write_timer(timers.add_and_start("write_output.write"));
  write_sdpb_input_files(output_dir, rank, num_procs, indices, objective_const,
                         dual_objectives_b, dual_constraint_groups);
  write_timer.stop();
}
