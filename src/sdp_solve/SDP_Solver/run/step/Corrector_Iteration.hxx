#pragma once

#include "sdpb_util/Boost_Float.hxx"

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
  // How fast mu decreases with time:
  // -d(log10 mu)/dt,
  //   where dt is time since solver step start.
  // Time t is normalized by full solver step (with single corrector iteration) timing.
  // In other words, log_mu_speed = 1
  //   if mu decreased by a factor of 10 after the first corrector iteration.
  Boost_Float log_mu_speed_full;
  // Same as log_mu_speed_full, but dt is current corrector iteration time.
  Boost_Float log_mu_speed_corrector;
  bool is_canceled = false;
};
