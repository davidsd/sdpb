#include "../Functional.hxx"

std::vector<El::BigFloat>
Functional::eval(const std::vector<El::BigFloat> &coords,
                 const std::vector<std::vector<El::BigFloat>> &optimals)
{
  std::vector<El::BigFloat> result;
  if(coords.size() != blocks.size() || coords.size() != optimals.size())
    {
      throw std::runtime_error("INTERNAL ERROR mismatch");
    }
  result.reserve(coords.size());
  for(size_t index(0); index != coords.size(); ++index)
    {
      if(has_prefactor)
        {
          throw std::runtime_error("unimpemented prefactor");
          // std::stringstream ss;
          // ss << coords[index];
          // Boost_Float x_mpfr(ss.str());
          // prefactor = El::BigFloat(
          //   to_string(pow(4 * (3 - sqrt(Boost_Float(2.0))), x_mpfr)));
        }
      result.emplace_back(blocks[index].eval(coords[index], optimals[index]));
    }
  return result;
}
