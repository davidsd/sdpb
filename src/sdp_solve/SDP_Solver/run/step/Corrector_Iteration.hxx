#pragma once

#include <El.hpp>

struct Corrector_Iteration final
{
  El::BigFloat primal_step_length, dual_step_length;
  El::BigFloat max_primal_step_length, max_dual_step_length;
  // mu' = Tr(X'Y')/X.dim
  // where
  //   X' = X + primal_step_length * dX
  //   Y' = Y + primal_step_length * dY
  El::BigFloat mu;
  // R_error = max(abs(mu'*I - X'Y'))
  El::BigFloat R_error;
  // R_mean_abs = mean(abs(mu'*I - X'Y'))
  El::BigFloat R_mean_abs;
  // 1 - min_step_length * (1 - beta)
  El::BigFloat reduce_factor;
};
