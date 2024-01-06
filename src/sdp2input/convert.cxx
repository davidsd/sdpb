#include "sdp_read/sdp_read.hxx"
#include "sdpb_util/Timers/Timers.hxx"
#include "sdpb_util/max_normalization_index.hxx"
#include "sdp_convert/sdp_convert.hxx"
#include "sdp_read/polynomial_matrices_with_prefactor/PMWP_SDP.hxx"

#include <filesystem>
#include <algorithm>

namespace fs = std::filesystem;

std::vector<Polynomial> bilinear_basis(const Damped_Rational &damped_rational,
                                       const size_t &half_max_degree);

namespace
{
  Polynomial_Vector_Matrix PMWP_to_PVM(
    const Positive_Matrix_With_Prefactor &positive_matrix_with_prefactor,
    const std::vector<El::BigFloat> &normalization,
    const size_t max_normalization_index, Timers &timers)
  {
    Polynomial_Vector_Matrix pvm;
    Scoped_Timer scalings_timer(timers, "scalings");
    const size_t max_degree([&]() {
      int64_t result(0);
      for(auto &pvv : positive_matrix_with_prefactor.polynomials)
        for(auto &pv : pvv)
          for(auto &polynomial : pv)
            {
              result = std::max(result, polynomial.degree());
            }
      return result;
    }());
    std::vector<Boost_Float> points(sample_points(max_degree + 1)),
      scalings(sample_scalings(
        points, positive_matrix_with_prefactor.damped_rational));
    scalings_timer.stop();

    pvm.rows = positive_matrix_with_prefactor.polynomials.size();
    pvm.cols = positive_matrix_with_prefactor.polynomials.front().size();

    Scoped_Timer bilinear_basis_timer(timers, "bilinear_basis");

    pvm.bilinear_basis = bilinear_basis(
      positive_matrix_with_prefactor.damped_rational, max_degree / 2);

    bilinear_basis_timer.stop();

    pvm.sample_points.reserve(points.size());
    for(auto &point : points)
      {
        pvm.sample_points.emplace_back(to_string(point));
      }
    pvm.sample_scalings.reserve(scalings.size());
    for(auto &scaling : scalings)
      {
        pvm.sample_scalings.emplace_back(to_string(scaling));
      }

    Scoped_Timer pvm_timer(timers, "pvm");
    pvm.elements.reserve(pvm.rows * pvm.cols);
    for(auto &pvv : positive_matrix_with_prefactor.polynomials)
      for(auto &pv : pvv)
        {
          pvm.elements.emplace_back();
          auto &pvm_polynomials(pvm.elements.back());
          pvm_polynomials.reserve(pv.size());
          pvm_polynomials.push_back(
            pv.at(max_normalization_index)
            / normalization.at(max_normalization_index));
          auto &pvm_constant(pvm_polynomials.back());

          for(size_t index = 0; index < normalization.size(); ++index)
            {
              if(index != max_normalization_index)
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
    return pvm;
  }
}

void convert(const PMWP_SDP &pwmp_sdp, El::BigFloat &objective_const,
             std::vector<El::BigFloat> &dual_objective_b,
             std::vector<Dual_Constraint_Group> &dual_constraint_groups,
             Timers &timers)
{
  Scoped_Timer convert_timer(timers, "convert");
  size_t max_index;
  {
    Scoped_Timer objectives_timer(timers, "objectives");
    max_index = max_normalization_index(pwmp_sdp.normalization);

    objective_const = pwmp_sdp.objective.at(max_index)
                      / pwmp_sdp.normalization.at(max_index);
    dual_objective_b = std::vector<El::BigFloat>();
    dual_objective_b.reserve(pwmp_sdp.normalization.size() - 1);
    for(size_t index = 0; index < pwmp_sdp.normalization.size(); ++index)
      {
        if(index != max_index)
          {
            dual_objective_b.push_back(pwmp_sdp.objective.at(index)
                                       - pwmp_sdp.normalization.at(index)
                                           * objective_const);
          }
      }
  }

  {
    Scoped_Timer matrices_timer(timers, "matrices");
    for(size_t i = 0; i < pwmp_sdp.matrices.size(); ++i)
      {
        size_t block_index = pwmp_sdp.matrix_index_local_to_global.at(i);
        const auto &pmwp_matrix = pwmp_sdp.matrices.at(i);
        Scoped_Timer index_timer(timers,
                                 "block_" + std::to_string(block_index));
        Polynomial_Vector_Matrix pvm = PMWP_to_PVM(
          pmwp_matrix, pwmp_sdp.normalization, max_index, timers);
        Scoped_Timer dual_constraint_timer(timers, "dual_constraint");
        dual_constraint_groups.emplace_back(block_index, std::move(pvm));
      }
  }
}
