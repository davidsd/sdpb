#include "sdp_solve/SDP_Solver.hxx"

// Step length is calculated as a function (see below)
// of max_step = \alpha(M, dM),
// which denotes the largest positive real number
// such that M + \alpha dM is positive semidefinite.
//
// \alpha(M, dM) is computed with a Cholesky decomposition M = L L^T.
// The eigenvalues of M + \alpha dM are equal to the eigenvalues of 1
// + \alpha L^{-1} dM L^{-T}.  The correct \alpha is then -1/lambda,
// where lambda is the smallest eigenvalue of L^{-1} dM L^{-T}.
//
// Inputs:
// - MCholesky = L, the Cholesky decomposition of M (M itself is not needed)
// - dM, a Block_Diagonal_Matrix with the same structure as M
// Workspace:
// - MInvDM (NB: overwritten when computing minEigenvalue)
// - eigenvalues, a Vector of eigenvalues for each block of M
// Output:
// - (output parameter:) max_step
// - (returned:) piecewise linear function of max_step = \alpha(M, dM),
// going through the points (0,0), (boost_step_min, gamma * boost_step_min), (boost_step_max, 1)
// and limited by 0 <= step <= 1.

// A := L^{-1} A L^{-T}
void lower_triangular_inverse_congruence(const Block_Diagonal_Matrix &L,
                                         Block_Diagonal_Matrix &A);

El::BigFloat min_eigenvalue(Block_Diagonal_Matrix &A);

// Compute step length for known max_step
El::BigFloat
step_length(const El::BigFloat &max_step, const El::BigFloat &gamma,
            const El::BigFloat &boost_step_min,
            const El::BigFloat &boost_step_max)
{
  if(max_step <= 0)
    {
      if(El::mpi::Rank() == 0)
        PRINT_WARNING(DEBUG_STRING(max_step));
      // TODO: in old algorithm, we return 1.
      // Shall we return 0 instead? Some tests fail (have slightly different output) if we do so.
      return 1;
    }

  // Old algorithm (disable boosting for large steps)
  if(boost_step_max < 0)
    return El::Min(gamma * max_step, El::BigFloat(1));

  ASSERT(boost_step_max >= boost_step_min, DEBUG_STRING(boost_step_min),
         DEBUG_STRING(boost_step_max));
  ASSERT(boost_step_min >= 0, DEBUG_STRING(boost_step_min));

  El::BigFloat step;
  if(max_step < boost_step_min)
    {
      // y = gamma * x
      step = gamma * max_step;
    }
  else if(max_step < boost_step_max)
    {
      // Linear function from (boost_step_min, gamma * boost_step_min) to (boost_step_max, 1)
      const El::BigFloat x1 = boost_step_min;
      const El::BigFloat x2 = boost_step_max;
      const El::BigFloat y1 = El::Min(gamma * x1, El::BigFloat(1));
      const El::BigFloat y2 = 1;

      step = y1 + (y2 - y1) / (x2 - x1) * (max_step - x1);
    }
  else
    {
      step = El::Max(max_step, El::BigFloat(1));
    }
  return El::Min(step, El::BigFloat(1));
}

// Old algorithm, without boosting
El::BigFloat
step_length(const El::BigFloat &max_step, const El::BigFloat &gamma)
{
  return step_length(max_step, gamma, -1, -1);
}

El::BigFloat
step_length(const Block_Diagonal_Matrix &MCholesky,
            const Block_Diagonal_Matrix &dM, const El::BigFloat &gamma,
            const El::BigFloat &boost_step_min,
            const El::BigFloat &boost_step_max, const std::string &timer_name,
            El::BigFloat &max_step, Timers &timers)
{
  Scoped_Timer step_length_timer(timers, timer_name);
  // MInvDM = L^{-1} dM L^{-T}, where M = L L^T
  Block_Diagonal_Matrix MInvDM(dM);
  lower_triangular_inverse_congruence(MCholesky, MInvDM);
  Scoped_Timer lambda_timer(timers, "min_eigenvalue");
  const El::BigFloat lambda(min_eigenvalue(MInvDM));
  lambda_timer.stop();

  // maxstep = \alpha from the comments above
  max_step = -1 / lambda;

  return step_length(max_step, gamma, boost_step_min, boost_step_max);
}
