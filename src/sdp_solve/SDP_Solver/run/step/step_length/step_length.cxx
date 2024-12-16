#include "sdp_solve/SDP_Solver.hxx"

// min(gamma \alpha(M, dM), 1), where \alpha(M, dM) denotes the
// largest positive real number such that M + \alpha dM is positive
// semidefinite.
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
// - min(\gamma \alpha(M, dM), 1) (returned)

// A := L^{-1} A L^{-T}
void lower_triangular_inverse_congruence(const Block_Diagonal_Matrix &L,
                                         Block_Diagonal_Matrix &A);

El::BigFloat min_eigenvalue(Block_Diagonal_Matrix &A);

El::BigFloat
step_length(const Block_Diagonal_Matrix &MCholesky,
            const Block_Diagonal_Matrix &dM, const El::BigFloat &gamma,
            const std::string &timer_name, El::BigFloat &max_step,
            Timers &timers)
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

  // Old algorithm:
  // return El::Min(gamma * maxstep, El::BigFloat(1));

  // New algorithm:

  // TODO does maxstep <= 0 make any sense?
  if(max_step >= 1 || max_step <= 0)
    return 1;

  // TODO this is some magic (by Ning). Need to explain and add new SDPB parameter.
  //
  // Old algorithm (with default gamma = --stepLengthReduction = 0.7):
  //   step(maxstep=0.95) = gamma * 0.95 = 0.665
  //   step(maxstep=1.00) = gamma = 0.7
  //   step(maxstep=1/gamma) = 1
  // New algorithm:
  //   step(maxstep=0.95) = 0.75
  //   step(maxstep=1.00) = 1
  //   step(maxstep=1/gamma) = 1
  // i.e. we increased step length for maxstep in range (0.95, 1/gamma) = (0.95, 1.43)
  if(max_step > 0.95)
    return 1 - 5 * (1 - max_step);
  return gamma * max_step;
}
