#include "../Boost_Float.hxx"
#include "../sdp_convert.hxx"
#include "../sdp_convert/write_vector.hxx"

void write_functions(
  const boost::filesystem::path &output_path,
  const std::vector<El::BigFloat> &dual_objective_b,
  const std::vector<Polynomial_Vector_Matrix> &polynomial_vector_matrices)
{
  boost::filesystem::ofstream output_stream(output_path);
  set_stream_precision(output_stream);
  output_stream << "{\n  \"objective\":\n";
  write_vector(output_stream, dual_objective_b, "    ");
  output_stream << ",\n  \"normalization\":\n";
  std::vector<int64_t> normalization(dual_objective_b.size(), 0);
  normalization.at(0) = 1;
  write_vector(output_stream, normalization, "  ");
  output_stream << ",\n  \"functions\":\n  [\n";

  const Boost_Float pi(boost::math::constants::pi<Boost_Float>());
  const El::BigFloat zero(0);

  size_t block_number(0);
  for(auto block(polynomial_vector_matrices.begin());
      block != polynomial_vector_matrices.end(); ++block)
    {
      if(block != polynomial_vector_matrices.begin())
        {
          output_stream << ",\n";
        }
      output_stream << "    [\n";

      const size_t num_chebyshev_points([&block]() {
        size_t max(0);
        for(auto &matrix_element : block->elements)
          for(auto &poly : matrix_element)
            {
              max = std::max(max, poly.coefficients.size());
            }
        return max;
      }());
      const El::BigFloat max_delta([&block]() {
        El::BigFloat max(0);
        for(auto &point : block->sample_points)
          {
            max = std::max(max, point);
          }
        return 8 * max;
      }());
      std::vector<El::BigFloat> chebyshev_zeros(num_chebyshev_points, 0);
      for(size_t index(0); index != chebyshev_zeros.size(); ++index)
        {
          chebyshev_zeros[index]
            = 0.5 * max_delta
              * El::BigFloat(
                to_string(1
                          + cos((pi * (num_chebyshev_points - index - 0.5))
                                / num_chebyshev_points)));
        }

      std::vector<int64_t> max_degree(block->rows * block->cols, 0),
        min_degree(block->rows * block->cols,
                   std::numeric_limits<int64_t>::max());
      for(int64_t row(0); row != block->rows; ++row)
        {
          for(int64_t column(0); column != block->cols; ++column)
            {
              for(auto &poly : block->elt(row, column))
                {
                  // Sometimes, the highest degree coefficients are zero
                  for(size_t index(0); index != poly.coefficients.size();
                      ++index)
                    {
                      const int64_t degree(poly.coefficients.size() - 1
                                           - index);
                      if(poly.coefficients.at(degree) != zero)
                        {
                          max_degree.at(row + column * block->rows) = std::max(
                            max_degree.at(row + column * block->rows), degree);
                          min_degree.at(row + column * block->rows) = std::min(
                            min_degree.at(row + column * block->rows), degree);
                        }
                    }
                }
            }
        }
      // Fix max_degree so that we have the correct limiting
      // determinant
      switch(block->rows)
        {
        case 1: break;
          case 2: {
            const int64_t first(max_degree.at(0) + max_degree.at(3)),
              second(2 * max_degree.at(1));
            if(first > second)
              {
                max_degree.at(1) = 0;
                max_degree.at(2) = 0;
              }
            else if(first < second)
              {
                max_degree.at(0) = 0;
                max_degree.at(3) = 0;
              }
          }
          break;
        default:
          throw std::runtime_error(
            "Too large a dimension.  Only 1x1 and 2x2 supported: "
            + std::to_string(block->rows));
        }

      for(int64_t row(0); row != block->rows; ++row)
        {
          if(row != 0)
            {
              output_stream << ",\n";
            }
          output_stream << "      [\n";
          for(int64_t column(0); column != block->cols; ++column)
            {
              if(column != 0)
                {
                  output_stream << ",\n";
                }
              output_stream << "        [\n";
              size_t poly_number(0);
              for(auto poly(block->elt(row, column).begin());
                  poly != block->elt(row, column).end(); ++poly)
                {
                  if(poly != block->elt(row, column).begin())
                    {
                      output_stream << ",\n";
                    }
                  output_stream << "          {\n"
                                << "            \"max_delta\": \"" << max_delta
                                << "\",\n"
                                << "            \"infinity_value\": \"";
                  if(poly->degree()
                     < max_degree.at(row + column * block->rows))
                    {
                      output_stream << "0";
                    }
                  else
                    {
                      output_stream << poly->coefficients.at(
                        max_degree.at(row + column * block->rows));
                    }
                  output_stream << "\",\n"
                                << "            \"epsilon_value\": \"";
                  if(poly->degree()
                     < min_degree.at(row + column * block->rows))
                    {
                      output_stream << "0";
                    }
                  else
                    {
                      output_stream << poly->coefficients.at(
                        min_degree.at(row + column * block->rows));
                    }
                  output_stream << "\",\n"
                                << "            \"chebyshev_values\":\n"
                                << "            [\n";

                  for(auto point(chebyshev_zeros.begin());
                      point != chebyshev_zeros.end(); ++point)
                    {
                      if(point != chebyshev_zeros.begin())
                        {
                          output_stream << ",\n";
                        }
                      output_stream << "              \"" << (*poly)(*point)
                                    << "\"";
                    }
                  output_stream << "\n            ]"
                                << "\n          }";
                  ++poly_number;
                }
              output_stream << "\n        ]";
            }
          output_stream << "\n      ]";
        }
      output_stream << "\n    ]";
      ++block_number;
    }
  output_stream << "\n  ]\n}\n";
}
