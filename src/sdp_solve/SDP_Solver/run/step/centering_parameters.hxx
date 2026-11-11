#pragma once

#include "sdp_solve/Block_Matrix/Abstract_Block_Diagonal_Matrix.hxx"
#include "sdp_solve/Solver_Parameters.hxx"

// Centering parameter \beta_p for the predictor step
inline El::BigFloat predictor_centering_parameter(const Solver_Parameters &parameters,
                                           const bool is_primal_dual_feasible)
{
  return is_primal_dual_feasible ? El::BigFloat(0)
                                 : parameters.infeasible_centering_parameter;
}

// Centering parameter \beta_c for the corrector step
template<class Derived>
El::BigFloat corrector_centering_parameter(
  const Solver_Parameters &parameters, const Abstract_Block_Diagonal_Matrix<Derived> &X,
  const Abstract_Block_Diagonal_Matrix<Derived> &dX, const Abstract_Block_Diagonal_Matrix<Derived> &Y,
  const Abstract_Block_Diagonal_Matrix<Derived> &dY, const El::BigFloat &mu,
  const bool is_primal_dual_feasible, const size_t &total_psd_rows)
{
  El::BigFloat r
    = frobenius_product_of_sums(X, dX, Y, dY) / (mu * total_psd_rows);
  El::BigFloat beta = r < 1 ? r * r : r;

  if(is_primal_dual_feasible)
    {
      return Min(Max(parameters.feasible_centering_parameter, beta),
                 El::BigFloat(1));
    }
  else
    {
      return Max(parameters.infeasible_centering_parameter, beta);
    }
}
