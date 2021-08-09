#include "../sdp_read.hxx"
#include "../sdp_convert/write_vector.hxx"

#include <boost/filesystem.hpp>

void write_functions(
  const boost::filesystem::path &output_path,
  const std::vector<El::BigFloat> &objectives,
  const std::vector<El::BigFloat> &normalization,
  const std::vector<Positive_Matrix_With_Prefactor> &matrices)
{
  boost::filesystem::ofstream output_stream(output_path);
  set_stream_precision(output_stream);
  output_stream << "{\n  \"objective\":\n";
  write_vector(output_stream, objectives, "    ");
  output_stream << ",\n  \"normalization\":\n";
  write_vector(output_stream, normalization, "  ");
  output_stream << ",\n  \"functions\":\n  [\n";

  const Boost_Float pi(boost::math::constants::pi<Boost_Float>());
  const El::BigFloat zero(0);

  size_t block_number(0);
  for(auto block(matrices.begin()); block != matrices.end(); ++block)
    {
      if(block != matrices.begin())
        {
          output_stream << ",\n";
        }
      output_stream << "    [\n";

      const size_t num_chebyshev_points([&block]() {
        size_t max(0);
        for(auto &row : block->polynomials)
          for(auto &column : row)
            for(auto &poly : column)
              {
                max = std::max(max, poly.coefficients.size());
              }
        return max;
      }());

      const El::BigFloat max_delta(6 * num_chebyshev_points);
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

      const size_t num_rows(block->polynomials.size()),
        num_columns(block->polynomials.front().size());
      std::vector<int64_t> max_degree(num_rows * num_columns, 0),
        min_degree(max_degree.size(), std::numeric_limits<int64_t>::max());
      for(size_t row(0); row != num_rows; ++row)
        {
          for(size_t column(0); column != num_columns; ++column)
            {
              for(auto &poly : block->polynomials[row][column])
                {
                  // Sometimes, the highest degree coefficients are zero
                  for(size_t index(0); index != poly.coefficients.size();
                      ++index)
                    {
                      const int64_t degree(poly.coefficients.size() - 1
                                           - index);
                      if(poly.coefficients.at(degree) != zero)
                        {
                          max_degree.at(row + column * num_rows) = std::max(
                            max_degree.at(row + column * num_rows), degree);
                          min_degree.at(row + column * num_rows) = std::min(
                            min_degree.at(row + column * num_rows), degree);
                        }
                    }
                }
            }
        }
      // Fix max_degree so that we have the correct limiting
      // determinant
      switch(num_rows)
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
            + std::to_string(num_rows));
        }

      for(size_t row(0); row != num_rows; ++row)
        {
          if(row != 0)
            {
              output_stream << ",\n";
            }
          output_stream << "      [\n";
          for(size_t column(0); column != num_columns; ++column)
            {
              if(column != 0)
                {
                  output_stream << ",\n";
                }
              output_stream << "        [\n";
              size_t poly_number(0);
              for(auto poly(block->polynomials[row][column].begin());
                  poly != block->polynomials[row][column].end(); ++poly)
                {
                  if(poly != block->polynomials[row][column].begin())
                    {
                      output_stream << ",\n";
                    }
                  output_stream << "          {\n"
                                << "            \"max_delta\": \"" << max_delta
                                << "\",\n"
                                << "            \"infinity_value\": \"";
                  if(poly->degree() < max_degree.at(row + column * num_rows))
                    {
                      output_stream << "0";
                    }
                  else
                    {
                      output_stream << poly->coefficients.at(
                        max_degree.at(row + column * num_rows));
                    }
                  output_stream << "\",\n"
                                << "            \"epsilon_value\": \"";
                  if(poly->degree() < min_degree.at(row + column * num_rows))
                    {
                      output_stream << "0";
                    }
                  else
                    {
                      output_stream << poly->coefficients.at(
                        min_degree.at(row + column * num_rows));
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
