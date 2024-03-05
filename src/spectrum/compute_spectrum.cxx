#include "eval_summed.hxx"
#include "get_zeros.hxx"
#include "Zeros.hxx"
#include "compute_lambda.hxx"
#include "sdpb_util/Mesh.hxx"
#include "pmp/max_normalization_index.hxx"
#include "sdpb_util/fill_weights.hxx"

std::vector<Zeros>
compute_spectrum(const Polynomial_Matrix_Program &pmp,
                 const El::Matrix<El::BigFloat> &y,
                 const std::vector<El::Matrix<El::BigFloat>> &x,
                 const El::BigFloat &threshold,
                 const El::BigFloat &mesh_threshold, const bool &need_lambda)
{
  std::vector<El::BigFloat> normalization;
  if(pmp.normalization.has_value())
    normalization = pmp.normalization.value();
  if(pmp.normalization.has_value())
    {
      normalization = pmp.normalization.value();
    }
  else
    {
      normalization.resize(pmp.objective.size(), 0);
      normalization.at(0) = 1;
    }

  const auto &matrices = pmp.matrices;

  const size_t max_index(max_normalization_index(normalization));
  std::vector<El::BigFloat> weights(normalization.size());
  fill_weights(y, max_index, normalization, weights);

  const El::BigFloat zero(0);
  std::vector<Zeros> zeros_blocks(matrices.size());
  for(size_t block_index(0); block_index < matrices.size(); ++block_index)
    {
      zeros_blocks.at(block_index).block_path
        = pmp.block_paths.at(block_index);
      auto &block(matrices[block_index]);
      const size_t max_number_terms([&block]() {
        size_t max(0);
        for(int i = 0; i < block.polynomials.Height(); ++i)
          for(int j = 0; j < block.polynomials.Width(); ++j)
            for(const auto &poly : block.polynomials(i, j))
              {
                max = std::max(max, poly.coefficients.size());
              }
        return max;
      }());

      // The factor of 6 comes from the limiting scale for Laguerre
      // polynomial roots.
      const El::BigFloat max_delta(6 * max_number_terms);

      const size_t num_rows(block.polynomials.Height()),
        num_columns(block.polynomials.Width());

      std::vector<Polynomial_Vector> summed_polynomials(num_rows);
      El::BigFloat block_scale(0), product;
      for(size_t row(0); row != num_rows; ++row)
        {
          summed_polynomials[row].reserve(num_columns);
          for(size_t column(0); column != num_columns; ++column)
            {
              auto &summed(summed_polynomials[row].emplace_back());
              summed.coefficients.resize(max_number_terms, zero);
              for(size_t dual_index(0);
                  dual_index != block.polynomials(row, column).size();
                  ++dual_index)
                {
                  auto &poly(block.polynomials(row, column)[dual_index]);
                  for(size_t coeff(0); coeff != poly.coefficients.size();
                      ++coeff)
                    {
                      product = weights[dual_index] * poly.coefficients[coeff];
                      block_scale = std::max(block_scale, El::Abs(product));
                      summed.coefficients[coeff] += product;
                    }
                }
            }
        }
      const El::BigFloat block_epsilon(block_scale
                                       * El::limits::Epsilon<El::BigFloat>());

      Mesh mesh(
        zero, max_delta,
        [&](const El::BigFloat &point) {
          return eval_summed(summed_polynomials, point);
        },
        mesh_threshold, block_epsilon);

      auto &zeros(zeros_blocks.at(block_index).zeros);
      if(need_lambda)
        {
          compute_lambda(block.sample_points, block.sample_scalings,
                         block.polynomials.Height(), x.at(block_index),
                         get_zeros(mesh, threshold), zeros,
                         zeros_blocks.at(block_index).error);
        }
      else
        {
          for(auto &zero : get_zeros(mesh, threshold))
            {
              zeros.emplace_back(zero);
            }
        }
    }
  return zeros_blocks;
}
