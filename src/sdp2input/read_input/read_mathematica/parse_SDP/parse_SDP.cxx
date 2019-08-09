#include "parse_vector.hxx"
#include "parse_generic.hxx"
#include "../../Positive_Matrix_With_Prefactor_State.hxx"

#include <algorithm>
#include <iterator>
#include <string>

const char *parse_SDP(const char *begin, const char *end,
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

  std::vector<El::BigFloat> temp_vector;
  const char *end_objective(parse_vector(SDP_start, end, temp_vector));
  if(!temp_vector.empty())
    {
      std::swap(objectives, temp_vector);
      temp_vector.clear();
    }

  const char *comma(std::find(end_objective, end, ','));
  if(comma == end)
    {
      throw std::runtime_error("Missing comma after objective");
    }

  parse_vector(std::next(comma), end, temp_vector);
  if(!temp_vector.empty())
    {
      std::swap(normalization, temp_vector);
      temp_vector.clear();
    }

  comma = std::find(end_objective, end, ',');
  if(comma == end)
    {
      throw std::runtime_error("Missing comma after normalization");
    }

  // TODO: Create dual_constraint_groups directly.  After some
  // experiments, this does not seem to help with memory.

  // TODO: Write everything but the objectives directly.  This would
  // help significantly with memory, but it requires having the global
  // normalization available in the first file we read.
  std::vector<Positive_Matrix_With_Prefactor> temp_matrices;
  const char *end_matrices(
    parse_generic(std::next(comma), end, temp_matrices));
  {
    size_t offset(matrices.size());
    matrices.resize(matrices.size() + temp_matrices.size());
    for(size_t index = 0; index < temp_matrices.size(); ++index)
      {
        std::swap(matrices[offset + index], temp_matrices[index]);
      }
  }
  const auto close_bracket(std::find(end_matrices, end, ']'));
  if(close_bracket == end)
    {
      throw std::runtime_error("Missing ']' at end of SDP");
    }
  return std::next(close_bracket);
}
