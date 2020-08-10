#include "poles_prefactor.hxx"
#include "power_prefactor.hxx"
#include "../sdp_read.hxx"

El::BigFloat
eval_weighted(const Positive_Matrix_With_Prefactor &matrix,
              const El::BigFloat &x, const std::vector<El::BigFloat> &weights)
{
  if(matrix.polynomials.size() != 1 || matrix.polynomials.front().size() != 1)
    {
      throw std::runtime_error("Only handling 1x1 for now");
    }

  auto &polys(matrix.polynomials.front().front());
  if(weights.size() != polys.size())
    {
      throw std::runtime_error("INTERNAL ERROR mismatch: "
                               + std::to_string(weights.size()) + " "
                               + std::to_string(polys.size()));
    }

  El::BigFloat result(0);
  for(size_t index(0); index != weights.size(); ++index)
    {
      result += weights[index] * polys[index](x);
    }
  if(!matrix.damped_rational.is_constant())
    {
      result *= poles_prefactor(matrix.damped_rational.poles, x)
                * power_prefactor(matrix.damped_rational.base, x);
    }

  return result;
}
