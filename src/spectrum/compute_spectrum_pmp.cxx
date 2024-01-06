#include "Zeros.hxx"
#include "eval_summed.hxx"
#include "get_zeros.hxx"
#include "compute_lambda.hxx"
#include "sdp_read/sdp_read.hxx"
#include "sdp_convert/sdp_convert.hxx"
#include "sdpb_util/Mesh.hxx"
#include "sdpb_util/max_normalization_index.hxx"
#include "sdpb_util/fill_weights.hxx"

std::vector<Zeros> compute_spectrum_pmp(
  const std::vector<El::BigFloat> &normalization,
  const El::Matrix<El::BigFloat> &y,
  const std::vector<Positive_Matrix_With_Prefactor> &matrices,
  const std::vector<size_t> &block_indices,
  const std::vector<El::Matrix<El::BigFloat>> &x,
  const El::BigFloat &threshold, El::BigFloat &mesh_threshold,
  const bool &need_lambda)
{
  const size_t max_index(max_normalization_index(normalization));
  std::vector<El::BigFloat> weights(normalization.size());
  fill_weights(y, max_index, normalization, weights);

  const El::BigFloat zero(0);
  std::vector<Zeros> zeros_blocks(matrices.size());
  for(size_t local_index = 0; local_index < matrices.size(); ++local_index)
    {
      const auto &block_index = block_indices.at(local_index);
      const auto &block = matrices.at(local_index);
      const size_t max_number_terms([&block]() {
        size_t max(0);
        for(auto &row : block.polynomials)
          for(auto &column : row)
            for(auto &poly : column)
              {
                max = std::max(max, poly.coefficients.size());
              }
        return max;
      }());

      // The factor of 6 comes from the limiting scale for Laguerre
      // polynomial roots.
      const El::BigFloat max_delta(6 * max_number_terms);

      const size_t num_rows(block.polynomials.size()),
        num_columns(block.polynomials.front().size());

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
                  dual_index != block.polynomials[row][column].size();
                  ++dual_index)
                {
                  auto &poly(block.polynomials[row][column][dual_index]);
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

      // TODO: This should all happen in Boost_Float so that we do not
      // have to do all of these conversions.
      std::stringstream ss;
      set_stream_precision(ss);
      Mesh mesh(
        zero, max_delta,
        [&](const El::BigFloat &x) {
          ss.str("");
          ss << x;
          Boost_Float x_Boost_Float(ss.str());
          Boost_Float numerator(
            block.damped_rational.constant
            * pow(block.damped_rational.base, x_Boost_Float));
          Boost_Float denominator(1);
          for(auto &pole : block.damped_rational.poles)
            {
              denominator *= (x_Boost_Float - pole);
            }
          return El::BigFloat(to_string(numerator / denominator))
                 * eval_summed(summed_polynomials, x);
        },
        mesh_threshold, block_epsilon);

      auto &zeros(zeros_blocks.at(local_index).zeros);
      if(need_lambda)
        {
          std::vector<Boost_Float> points_boost(
            sample_points(max_number_terms)),
            scalings_boost(
              sample_scalings(points_boost, block.damped_rational));
          std::vector<El::BigFloat> points_gmp, scalings_gmp;
          for(auto &point : points_boost)
            {
              points_gmp.emplace_back(to_string(point));
            }
          for(auto &scaling : scalings_boost)
            {
              scalings_gmp.emplace_back(to_string(scaling));
            }
          compute_lambda(points_gmp, scalings_gmp, num_rows, x.at(local_index),
                         get_zeros(mesh, threshold), zeros,
                         zeros_blocks.at(local_index).error);
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
