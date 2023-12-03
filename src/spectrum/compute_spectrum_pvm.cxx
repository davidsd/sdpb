#include "eval_summed.hxx"
#include "get_zeros.hxx"
#include "Zeros.hxx"
#include "compute_lambda.hxx"
#include "sdpb_util/Mesh.hxx"
#include "sdpb_util/max_normalization_index.hxx"
#include "sdpb_util/fill_weights.hxx"

std::vector<Zeros>
compute_spectrum_pvm(const El::Matrix<El::BigFloat> &y,
                     const std::vector<Polynomial_Vector_Matrix> &matrices,
                     const std::vector<El::Matrix<El::BigFloat>> &x,
                     const El::BigFloat &threshold,
                     El::BigFloat &mesh_threshold, const bool &need_lambda)
{
  // pvm2sdp implicitly uses the first element as the normalized column
  std::vector<El::BigFloat> normalization(y.Height() + 1, 0);
  normalization.at(0) = 1;
  const size_t max_index(0);
  std::vector<El::BigFloat> weights(normalization.size());
  fill_weights(y, max_index, normalization, weights);

  const El::BigFloat zero(0);
  std::vector<Zeros> zeros_blocks(matrices.size());
  for(size_t block_index(0); block_index < matrices.size(); ++block_index)
    {
      auto &block(matrices[block_index]);
      const size_t max_number_terms([&block]() {
        size_t max(0);
        for(auto &element : block.elements)
          for(auto &poly : element)
            {
              max = std::max(max, poly.coefficients.size());
            }
        return max;
      }());

      // The factor of 6 comes from the limiting scale for Laguerre
      // polynomial roots.
      const El::BigFloat max_delta(6 * max_number_terms);

      const size_t num_rows(block.rows), num_columns(block.cols);

      std::vector<std::vector<Polynomial>> summed_polynomials(num_rows);
      El::BigFloat block_scale(0), product;
      for(size_t row(0); row != num_rows; ++row)
        {
          summed_polynomials[row].reserve(num_columns);
          for(size_t column(0); column != num_columns; ++column)
            {
              auto &summed(summed_polynomials[row].emplace_back());
              summed.coefficients.resize(max_number_terms, zero);
              for(size_t dual_index(0);
                  dual_index != block.elt(row, column).size(); ++dual_index)
                {
                  auto &poly(block.elt(row, column)[dual_index]);
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
                         block.rows, x.at(block_index),
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
