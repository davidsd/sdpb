#pragma once

#include "step_length_double.hxx"
#include "sdp_solve/Block_Matrix/Abstract_Block_Diagonal_Matrix.hxx"
#include "sdpb_util/Timers/Timers.hxx"

#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>
#include <sstream>

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

// Options of step_length(). Default: El::HermitianEig on every
// block.
struct Step_Length_Options
{
  // --stepLengthDouble: each block's smallest eigenvalue in double, and
  // El::HermitianEig only for a block whose double eigenvalue could move the
  // step length by more than double_tol (step_length_double.hxx; the GPU
  // backend's scheme)
  bool use_double = false;
  double double_tol = 1e-10; // --stepLengthDoubleTol
  // --stepLengthCheck: compute both, use El::HermitianEig's (the run is the
  // default one), print the agreement on rank 0
  bool check = false;
  // print the counts on rank 0 (--verbosity debug)
  bool print_stats = false;
};

// Smallest eigenvalue of one block by El::HermitianEig.
// Annoyingly, El::HermitianEig modifies 'block'.  It is OK, because
// it is only called from step_length(), which passes in a temporary.
inline El::BigFloat block_min_eigenvalue(El::DistMatrix<El::BigFloat> &block)
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
  hermitian_eig_ctrl.tridiagEigCtrl.dcCtrl.secularCtrl.maxIterations = 16384;
  El::HermitianEig(El::UpperOrLowerNS::LOWER, block, eigenvalues,
                   hermitian_eig_ctrl);
  return El::Min(eigenvalues);
}

// Minimum eigenvalue of A.  A is assumed to be symmetric.
template <class Derived>
El::BigFloat min_eigenvalue(Abstract_Block_Diagonal_Matrix<Derived> &A)
{
  El::BigFloat local_min(El::limits::Max<El::BigFloat>());

  for(auto &block : A.blocks)
    {
      local_min = El::Min(local_min, block_min_eigenvalue(block));
    }
  return El::mpi::AllReduce(local_min, El::mpi::MIN, El::mpi::COMM_WORLD);
}

inline El::BigFloat
step_from_eigenvalue(const El::BigFloat &lambda, const El::BigFloat &gamma)
{
  if(lambda > -gamma)
    return 1;
  return -gamma / lambda;
}

// log10 |x| (-inf for 0)
inline double step_length_log10_abs(const El::BigFloat &x)
{
  if(x.Sign() == 0)
    return -std::numeric_limits<double>::infinity();
  long ex;
  const double d = mpf_get_d_2exp(&ex, x.gmp_float.get_mpf_t());
  return std::log10(std::fabs(d)) + ex * std::log10(2.0);
}

// digits of agreement of a with the reference b: -log10 |a - b| / |b|
inline double step_length_digits(const El::BigFloat &a, const El::BigFloat &b)
{
  const El::BigFloat d = El::Abs(a - b);
  if(d.Sign() == 0)
    return std::numeric_limits<double>::infinity();
  return step_length_log10_abs(b) - step_length_log10_abs(d);
}

// --stepLengthDouble: the minimum over blocks of each block's smallest
// eigenvalue, in double where that resolves the block's step length to
// options.double_tol, by El::HermitianEig otherwise. Zero blocks are
// ignored (+Max). With options.check, returns El::HermitianEig's minimum (as
// min_eigenvalue()) and prints how the two agree.
template <class Derived>
El::BigFloat
min_eigenvalue_double(Abstract_Block_Diagonal_Matrix<Derived> &A,
                      const El::BigFloat &gamma,
                      const Step_Length_Options &options,
                      const std::string &name)
{
  const El::BigFloat max_value(El::limits::Max<El::BigFloat>());
  El::BigFloat local_min(max_value), local_stock_min(max_value);
  double n_blocks = 0, n_accepted = 0, n_flagged = 0, n_ignored = 0;
  // check: over the accepted blocks, the worst relative error of the block's
  // step length and the worst agreement of its eigenvalue (digits)
  El::BigFloat max_step_error(0);
  double min_lambda_digits = std::numeric_limits<double>::infinity();
  const double gamma_d = double(gamma);
  for(auto &block : A.blocks)
    {
      ++n_blocks;
      const Double_Min_Eig r
        = double_min_eigenvalue(block, gamma_d, options.double_tol);
      // every rank of the grid has the same double matrix and makes the same
      // decision; make sure, since El::HermitianEig is collective
      int flagged = r.status == Double_Min_Eig_Status::flagged ? 1 : 0;
      if(block.Grid().Size() > 1)
        flagged
          = El::mpi::AllReduce(flagged, El::mpi::MAX, block.DistComm());
      if(flagged)
        {
          ++n_flagged;
          const El::BigFloat value = block_min_eigenvalue(block);
          local_min = El::Min(local_min, value);
          local_stock_min = El::Min(local_stock_min, value);
          continue;
        }
      if(r.status == Double_Min_Eig_Status::accepted)
        {
          ++n_accepted;
          local_min = El::Min(local_min, r.value);
        }
      else
        ++n_ignored;
      if(options.check)
        {
          const El::BigFloat stock = block_min_eigenvalue(block);
          local_stock_min = El::Min(local_stock_min, stock);
          if(r.status == Double_Min_Eig_Status::accepted)
            {
              const El::BigFloat s_stock = step_from_eigenvalue(stock, gamma);
              const El::BigFloat s_double
                = step_from_eigenvalue(r.value, gamma);
              max_step_error = El::Max(max_step_error,
                                       El::Abs(s_double - s_stock) / s_stock);
              min_lambda_digits = std::min(
                min_lambda_digits, step_length_digits(r.value, stock));
            }
        }
    }
  const El::BigFloat global_min
    = El::mpi::AllReduce(local_min, El::mpi::MIN, El::mpi::COMM_WORLD);
  if(!(options.print_stats || options.check))
    return global_min;

  // only the first rank of each grid counts its blocks
  const bool counts
    = A.blocks.empty() || A.blocks.front().Grid().Rank() == 0;
  std::vector<double> sums
    = {counts ? n_blocks : 0., counts ? n_accepted : 0.,
       counts ? n_flagged : 0., counts ? n_ignored : 0.};
  El::mpi::AllReduce(sums.data(), sums.size(), El::mpi::SUM,
                     El::mpi::COMM_WORLD);
  std::ostringstream os;
  os << "stepLengthDouble " << name << ": blocks " << sums[0]
     << " double " << sums[1] << " multiprecision " << sums[2] << " zero "
     << sums[3];
  El::BigFloat stock_min;
  if(options.check)
    {
      stock_min = El::mpi::AllReduce(local_stock_min, El::mpi::MIN,
                                     El::mpi::COMM_WORLD);
      max_step_error = El::mpi::AllReduce(max_step_error, El::mpi::MAX,
                                          El::mpi::COMM_WORLD);
      min_lambda_digits = El::mpi::AllReduce(
        min_lambda_digits, El::mpi::MIN, El::mpi::COMM_WORLD);
      const El::BigFloat step_stock = step_from_eigenvalue(stock_min, gamma);
      const El::BigFloat step_double = step_from_eigenvalue(global_min, gamma);
      os << std::setprecision(4)
         << " | check: accepted blocks' max relative step error "
         << double(max_step_error) << ", min lambda digits "
         << min_lambda_digits << "; lambda_min " << std::setprecision(17)
         << double(stock_min) << " (default) step " << double(step_stock)
         << std::setprecision(4) << ", step relative error "
         << double(El::Abs(step_double - step_stock) / step_stock)
         << ", lambda digits " << step_length_digits(global_min, stock_min);
    }
  if(El::mpi::Rank() == 0)
    std::cout << os.str() << std::endl;
  return options.check ? stock_min : global_min;
}

template <class Derived>
El::BigFloat
step_length(const Abstract_Block_Diagonal_Matrix<Derived> &MCholesky,
            const Abstract_Block_Diagonal_Matrix<Derived> &dM,
            const El::BigFloat &gamma, const std::string &timer_name,
            Timers &timers,
            const Step_Length_Options &options = Step_Length_Options())
{
  Scoped_Timer step_length_timer(timers, timer_name);
  // MInvDM = L^{-1} dM L^{-T}, where M = L L^T
  if(!(options.use_double || options.check))
    {
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
  // --stepLengthDouble, --stepLengthCheck (sub-timers for the profile)
  Scoped_Timer congruence_timer(timers, "congruence");
  auto MInvDM(dM);
  lower_triangular_inverse_congruence(MCholesky, MInvDM);
  congruence_timer.stop();
  Scoped_Timer eig_timer(timers, "min_eigenvalue");
  const El::BigFloat lambda(
    min_eigenvalue_double(MInvDM, gamma, options, timer_name));
  eig_timer.stop();
  if(lambda > -gamma)
    {
      return 1;
    }
  else
    {
      return -gamma / lambda;
    }
}
