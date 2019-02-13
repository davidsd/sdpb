#include "parse_vector.hxx"
#include "../../Positive_Matrix_With_Prefactor_State.hxx"

#include <algorithm>
#include <iterator>
#include <string>

std::vector<char>::const_iterator
parse_SDP(const std::vector<char>::const_iterator &begin,
          const std::vector<char>::const_iterator &end,
          std::vector<El::BigFloat> &objectives,
          std::vector<El::BigFloat> &normalization,
          std::vector<Positive_Matrix_With_Prefactor> &matrices)

{
  const std::string SDP_literal("SDP[");
  auto SDP_start(
    std::search(begin, end, SDP_literal.begin(), SDP_literal.end()));
  if(SDP_start == end)
    {
      throw std::runtime_error("Could not find 'SDP['");
    }

  auto end_objective(parse_vector(SDP_start, end, objectives));

  auto comma(std::find(end_objective, end, ','));
  if(comma == end)
    {
      throw std::runtime_error("Missing comma after objective");
    }

  auto end_normalization(parse_vector(std::next(comma, 1), end, normalization));

  comma = std::find(end_objective, end, ',');
  if(comma == end)
    {
      throw std::runtime_error("Missing comma after normalization");
    }

  std::cout.precision(std::ceil(El::gmp::Precision() * std::log10(2.0)) + 1);
  for(auto &n: normalization)
    std::cout << "n: " << n << "\n";
  
  // auto end_matrices(parse_matrices(std::next(comma, 1), end));
  // const auto close_bracket(std::find(std::next(end_matrices, 1), end, ']'));
  // if(close_bracket == end)
  //   {
  //     throw std::runtime_error("Missing ']' at end of SDP");
  //   }
  // return std::next(close_bracket, 1);

  return end;
}
