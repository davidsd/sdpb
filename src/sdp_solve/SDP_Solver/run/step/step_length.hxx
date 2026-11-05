#pragma once

#include "sdp_solve/Block_Matrix/Abstract_Block_Diagonal_Matrix.hxx"
#include "sdpb_util/Timers/Timers.hxx"

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


// Minimum eigenvalue of A.  A is assumed to be symmetric.
// Annoyingly, El::HermitianEig modifies 'block'.  It is OK, because
// it is only called from step_length(), which passes in a temporary.
// Still ugly.
template <class Derived>
El::BigFloat min_eigenvalue(Abstract_Block_Diagonal_Matrix<Derived> &A)
{
  El::BigFloat local_min(El::limits::Max<El::BigFloat>());

  for(auto &block : A.blocks)
    {
      El::DistMatrix<El::BigFloat, El::VR, El::STAR> eigenvalues(block.Grid());
      /// There is a bug in El::HermitianEig when there is more than
      /// one level of recursion when computing eigenvalues.  One fix
      /// is to increase the cutoff so that there is no more than one
      /// level of recursion.

      /// An alternate workaround is to compute both eigenvalues and
      /// eigenvectors, but that seems to be significantly slower.
      El::HermitianEigCtrl<El::BigFloat> hermitian_eig_ctrl;
      hermitian_eig_ctrl.tridiagEigCtrl.dcCtrl.cutoff = block.Height() / 2 + 1;

      /// The default number of iterations is 40.  That is sometimes
      /// not enough, so we bump it up significantly.
      hermitian_eig_ctrl.tridiagEigCtrl.dcCtrl.secularCtrl.maxIterations
        = 16384;
      El::HermitianEig(El::UpperOrLowerNS::LOWER, block, eigenvalues,
                       hermitian_eig_ctrl);
      local_min = El::Min(local_min, El::Min(eigenvalues));
    }
  return El::mpi::AllReduce(local_min, El::mpi::MIN, El::mpi::COMM_WORLD);
}

template<class Derived>
El::BigFloat
step_length(const Abstract_Block_Diagonal_Matrix<Derived> &MCholesky,
            const Abstract_Block_Diagonal_Matrix<Derived> &dM, const El::BigFloat &gamma,
            const std::string &timer_name, Timers &timers)
{
  Scoped_Timer step_length_timer(timers, timer_name);
  // MInvDM = L^{-1} dM L^{-T}, where M = L L^T
  auto MInvDM(dM);
  lower_triangular_inverse_congruence(MCholesky, MInvDM);
  const El::BigFloat lambda(min_eigenvalue(MInvDM));
  if(lambda > -gamma)
    {
      return 1;
    }
  else
    {
      return -gamma / lambda;
    }
}
