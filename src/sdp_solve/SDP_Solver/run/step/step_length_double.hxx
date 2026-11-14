#pragma once

#include <El.hpp>

// --stepLengthDouble: the CPU port of the GPU backend's step eigenvalue
// (Gpu_Solver::min_eig, with Gpu_Compute's multiprecision fallback and
// its tolerance SDPB_GPU_EIG_TOL).
//
// For one block A of L^{-1} dM L^{-T}: A 2^{-e} in double (e = the largest
// exponent of the block's entries), symmetrized as (A + A^T)/2, its
// eigenvalues by LAPACK, lambda = the smallest, with the error bound
// err = (n + 2) 2^{-52} ||A||_F. The block's step length is
// step(lambda) = 1 if lambda > -gamma, else -gamma/lambda. If
// step(lambda - err) >= (1 - tol) step(lambda + err), the double lambda
// cannot move the step length by more than tol and is used; otherwise the
// block is flagged, and the caller computes its smallest eigenvalue with
// El::HermitianEig, as the default does. A zero block is ignored.

enum class Double_Min_Eig_Status
{
  ignored,  // empty or zero block
  accepted, // value = the double lambda (exactly, as a BigFloat)
  flagged   // the caller must use El::HermitianEig
};

struct Double_Min_Eig
{
  Double_Min_Eig_Status status = Double_Min_Eig_Status::ignored;
  El::BigFloat value;
};

// Collective over the block's grid. block is not changed.
Double_Min_Eig double_min_eigenvalue(const El::DistMatrix<El::BigFloat> &block,
                                     double gamma, double tol);
