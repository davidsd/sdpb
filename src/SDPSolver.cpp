//=======================================================================
// Copyright 2014-2015 David Simmons-Duffin.
// Distributed under the MIT License.
// (See accompanying file LICENSE or copy at
//  http://opensource.org/licenses/MIT)
//=======================================================================


#include <algorithm>
#include <iostream>
#include <vector>
#include "omp.h"
//Tweak to allow Ubuntu-14.04/gcc-4.8.4 and similar environments to compile
#define BOOST_NO_CXX11_SCOPED_ENUMS
#include <boost/filesystem.hpp>
#undef BOOST_NO_CXX11_SCOPED_ENUMS
#include "SDPSolver.h"
#include "Timers.h"

using boost::filesystem::path;
using boost::timer::nanosecond_type;
using std::cout;

// A note on function conventions: For functions which return large
// objects, we usually use references to the returned objects as
// arguments to a void function, which are then modified in place:
//
// void myFunction(const InputType1 &input1,
//                 const InputType2 &input2,
//                 ...
//                 OutputType1 &output1,
//                 OutputType2 &output2,
//                 ...) {
//   use input1, etc. to modify output1, etc.
// }
//
// We try to mark large inputs as `const' references.  However,
// several MBLAS and MLAPACK functions are not correctly annotated
// with `const's, so we occasionally have input arguments which cannot
// be marked as const because they're used in MBLAS/MLAPACK functions.
// We try to use comments to distinguish which arguments should be
// considered inputs and which should be considered outputs.

/***********************************************************************/
// Create and initialize an SDPSolver for the given SDP and
// SDPSolverParameters

SDPSolver::SDPSolver(const SDP &sdp, const SDPSolverParameters &parameters):
  sdp(sdp),
  parameters(parameters),
  x(sdp.primalObjective.size(), 0),
  X(sdp.psdMatrixBlockDims()),
  y(sdp.dualObjective.size(), 0),
  Y(X),
  dx(x),
  dX(X),
  dy(y),
  dY(Y),
  PrimalResidues(X),
  dualResidues(x),
  XCholesky(X),
  YCholesky(X),
  Z(X),
  R(X),
  BilinearPairingsXInv(sdp.bilinearPairingBlockDims()),
  BilinearPairingsY(BilinearPairingsXInv),
  SchurComplement(sdp.schurBlockDims()),
  SchurComplementCholesky(SchurComplement),
  SchurOffDiagonal(sdp.FreeVarMatrix),
  schurStabilizeIndices(SchurComplement.blocks.size()),
  schurStabilizeLambdas(SchurComplement.blocks.size()),
  stabilizeBlocks(SchurComplement.blocks.size()),
  Q(sdp.FreeVarMatrix.cols, sdp.FreeVarMatrix.cols),
  Qpivots(sdp.FreeVarMatrix.cols),
  dyExtended(Q.rows),
  StepMatrixWorkspace(X)
{
  // initialize bilinearPairingsWorkspace, eigenvaluesWorkspace, QRWorkspace
  for (unsigned int b = 0; b < sdp.bilinearBases.size(); b++) {
    bilinearPairingsWorkspace.push_back(Matrix(X.blocks[b].rows,
                                               BilinearPairingsXInv.blocks[b].cols));
    eigenvaluesWorkspace.push_back(Vector(X.blocks[b].rows));
    QRWorkspace.push_back(Vector(3*X.blocks[b].rows - 1));
  }

  // X = \Omega_p I
  X.addDiagonal(parameters.initialMatrixScalePrimal);
  // Y = \Omega_d I
  Y.addDiagonal(parameters.initialMatrixScaleDual);
}

/***********************************************************************/
// Specialized linear algebra helper functions needed for the solver.


// Result = Q'^T A Q', where Q' = Q \otimes 1, where \otimes denotes
// tensor product and 1 is an mxm identity matrix.
// Inputs:
// - A      : l*m x l*m symmetric matrix
// - Q      : l   x n   matrix
// Workspace:
// - Work   : l*m x n*m matrix, intermediate workspace (overwritten)
// Output:
// - Result : n*m x n*m symmetric matrix (overwritten)
//
// An explanation of the name: a 'congruence' refers to the action M
// -> A M A^T.  We use 'transpose congruence' to refer to a congruence
// with the transposed matrix M -> A^T M A.  'tensor' here refers to
// the fact that we're performing a congruence with the tensor product
// Q \otimes 1.
//
void tensorTransposeCongruence(const Matrix &A,
                               const Matrix &Q,
                               Matrix &Work,
                               Matrix &Result) {
  int m = A.rows / Q.rows;

  assert(Result.rows == Q.cols * m);
  assert(Result.cols == Q.cols * m);

  assert(Work.rows == A.rows);
  assert(Work.cols == Result.cols);

  // Work = A Q'
  for (int c = 0; c < Work.cols; c++) {
    int qCol       = c % Q.cols;
    int aColOffset = (c / Q.cols) * Q.rows;

    for (int r = 0; r < Work.rows; r++) {
      Real tmp = 0;
      for (int k = 0; k < Q.rows; k++) {
        tmp += A.elt(r, aColOffset + k) * Q.elt(k, qCol);
      }

      Work.elt(r, c) = tmp;
    }
  }

  // Result = Q'^T Work
  for (int c = 0; c < Result.cols; c++) {
    // since Result is symmetric, only compute its upper triangle
    for (int r = 0; r <= c; r++) {
      int qCol          = r % Q.cols;
      int workRowOffset = (r / Q.cols) * Q.rows;

      Real tmp = 0;
      for (int k = 0; k < Q.rows; k++) {
        tmp += Q.elt(k, qCol) * Work.elt(workRowOffset + k, c);
      }

      Result.elt(r, c) = tmp;

      // lower triangle is the same as upper triangle
      if (c != r) {
        Result.elt(c, r) = tmp;
      }
    }
  }
}

// Result_b = Q[b]'^T A_b Q[b]' for each block 0 <= b < Q.size()
// - Result_b, A_b denote the b-th blocks of Result, A, resp.
// - Q[b]' = Q[b] \otimes 1, where \otimes denotes tensor product
// - for each b, L.blocks[b], Q[b], Work[b], and Result.blocks[b] must
//   have the structure described above for
//   `tensorTransposeCongruence'
//
void blockTensorTransposeCongruence(const BlockDiagonalMatrix &A,
                                    const vector<Matrix> &Q,
                                    vector<Matrix> &Work,
                                    BlockDiagonalMatrix &Result) {
  #pragma omp parallel for schedule(dynamic)
  for (unsigned int b = 0; b < Q.size(); b++)
    tensorTransposeCongruence(A.blocks[b], Q[b], Work[b], Result.blocks[b]);
}

// Result = Q'^T A^{-1} Q', where Q' = Q \otimes 1, where \otimes
// denotes tensor product and 1 is an mxm idenity matrix.
// Inputs:
// - L      : l*m x l*m cholesky decomposition of A
// - Q      : l   x n   matrix
// Workspace:
// - Work   : l*m x n*m matrix, intermediate workspace (overwritten)
// Output:
// - Result : n*m x n*m symmetric matrix (overwritten)
//
void tensorInvTransposeCongruenceWithCholesky(const Matrix &L,
                                              const Matrix &Q,
                                              Matrix &Work,
                                              Matrix &Result) {
  // Work = L^{-1} (Q \otimes 1);
  for (int cw = 0; cw < Work.cols; cw++) {
    int mc  = cw / Q.cols;

    for (int rw = mc*Q.rows; rw < Work.rows; rw++) {
      int mr = rw / Q.rows;

      Real tmp = (mr != mc) ? Real(0) : Q.elt(rw % Q.rows, cw % Q.cols);
      for (int cl = mc*Q.rows; cl < rw; cl++)
        tmp -= L.elt(rw, cl)*Work.elt(cl, cw);

      Work.elt(rw, cw) = tmp/L.elt(rw, rw);
    }
  }

  // Result = Work^T Work
  for (int cr = 0; cr < Result.cols; cr++) {
    int mc = cr / Q.cols;

    for (int rr = 0; rr <= cr; rr++) {
      int mr = rr / Q.cols;

      Real tmp = 0;
      for (int rw = max(mr, mc)*Q.rows; rw < Work.rows; rw++)
        tmp += Work.elt(rw, cr)*Work.elt(rw, rr);

      Result.elt(rr, cr) = tmp;
      if (rr != cr)
        Result.elt(cr, rr) = tmp;
    }
  }
}

// Result_b = Q[b]'^T A_b^{-1} Q[b]' for each block 0 <= b < Q.size()
// - Result_b, A_b denote the b-th blocks of Result, A, resp.
// - Q[b]' = Q[b] \otimes 1, where \otimes denotes tensor product
// - for each b, L.blocks[b], Q[b], Work[b], and Result.blocks[b] must
//   have the structure described above for
//   `tensorInvTransposeCongruenceWithCholesky'
// 
void blockTensorInvTransposeCongruenceWithCholesky(const BlockDiagonalMatrix &L,
                                                   const vector<Matrix> &Q,
                                                   vector<Matrix> &Work,
                                                   BlockDiagonalMatrix &Result) {
  #pragma omp parallel for schedule(dynamic)
  for (unsigned int b = 0; b < Q.size(); b++)
    tensorInvTransposeCongruenceWithCholesky(L.blocks[b], Q[b],
                                             Work[b], Result.blocks[b]);
}

// Result^(blockRow,blockCol) = V D V^T, where D=diag(d) is a diagonal
// matrix.
//
// Here, we view Result as a matrix on the tensor product
// R^V.rows \otimes R^k.  Result^(blockRow,blockCol) refers to the
// block submatrix of size V.rows x V.rows at position (blockRow,
// blockCol) in the second tensor factor.
//
// Inputs:
// - d        : pointer to beginning of a length-V.cols vector
// - V        : V.rows x V.cols Matrix
// - blockRow : integer < k
// - blockCol : integer < k
// Output:
// - Result   : (k*V.rows) x (k*V.rows) square Matrix (overwritten)
//
void diagonalCongruence(Real const *d,
                        const Matrix &V,
                        const int blockRow,
                        const int blockCol,
                        Matrix &Result) {
  for (int p = 0; p < V.rows; p++) {
    for (int q = 0; q <= p; q++) {
      Real tmp = 0;

      for (int n = 0; n < V.cols; n++)
        tmp += *(d+n) * V.elt(p, n)*V.elt(q, n);

      Result.elt(blockRow*V.rows + p, blockCol*V.rows + q) = tmp;
      if (p != q)
        Result.elt(blockRow*V.rows + q, blockCol*V.rows + p) = tmp;
    }
  }
}

// v^T A^(blockRow, blockCol) v, where A^(r,s) is the (r,s)-th dim x
// dim block inside A.
//
// Input:
// - v        : pointer to the beginning of a vector of length dim
// - dim      : length of the vector v
// - A        : (k*dim) x (k*dim) matrix, where k > blockRow, blockCol
// - blockRow : integer labeling block of A
// - blockCol : integer labeling block of A
// Output: v^T A^(blockRow, blockCol) v (returned)
//
Real bilinearBlockPairing(const Real *v,
                          const int dim,
                          const Matrix &A,
                          const int blockRow,
                          const int blockCol) {
  Real result = 0;

  for (int r = 0; r < dim; r++) {
    Real tmp = 0;

    for (int c = 0; c < dim; c++)
      tmp += *(v+c) * A.elt(blockRow*dim + r, blockCol*dim + c);
    result += *(v+r) * tmp;
  }
  return result;
}

/***********************************************************************/
// Subroutines needed for each solver iteration

// Compute the SchurComplement matrix using BilinearPairingsXInv and
// BilinearPairingsY and the formula
//
//   S_{(j,r1,s1,k1), (j,r2,s2,k2)} = \sum_{b \in blocks[j]}
//          (1/4) (BilinearPairingsXInv_{ej s1 + k1, ej r2 + k2}*
//                 BilinearPairingsY_{ej s2 + k2, ej r1 + k1} +
//                 swaps (r1 <-> s1) and (r2 <-> s2))
// 
// where ej = d_j + 1.
//
// Inputs: sdp, BilinearPairingsXInv, BilinearPairingsY
// Output: SchurComplement (overwritten) 
//
void computeSchurComplement(const SDP &sdp,
                            const BlockDiagonalMatrix &BilinearPairingsXInv,
                            const BlockDiagonalMatrix &BilinearPairingsY,
                            BlockDiagonalMatrix &SchurComplement) {
  #pragma omp parallel for schedule(dynamic)
  for (unsigned int j = 0; j < sdp.dimensions.size(); j++) {
    const int ej = sdp.degrees[j] + 1;

    for (unsigned int u1 = 0; u1 < sdp.constraintIndices[j].size(); u1++) {
      const int ej_r1 = sdp.constraintIndices[j][u1].r * ej;
      const int ej_s1 = sdp.constraintIndices[j][u1].s * ej;
      const int k1    = sdp.constraintIndices[j][u1].k;

      for (unsigned int u2 = 0; u2 <= u1; u2++) {
        const int ej_r2 = sdp.constraintIndices[j][u2].r * ej;
        const int ej_s2 = sdp.constraintIndices[j][u2].s * ej;
        const int k2    = sdp.constraintIndices[j][u2].k;

        Real tmp = 0;
        for (vector<int>::const_iterator b = sdp.blocks[j].begin();
             b != sdp.blocks[j].end(); b++) {
          tmp += (BilinearPairingsXInv.blocks[*b].elt(ej_s1 + k1, ej_r2 + k2) *
                  BilinearPairingsY   .blocks[*b].elt(ej_s2 + k2, ej_r1 + k1) +
                  BilinearPairingsXInv.blocks[*b].elt(ej_r1 + k1, ej_r2 + k2) *
                  BilinearPairingsY   .blocks[*b].elt(ej_s2 + k2, ej_s1 + k1) +
                  BilinearPairingsXInv.blocks[*b].elt(ej_s1 + k1, ej_s2 + k2) *
                  BilinearPairingsY   .blocks[*b].elt(ej_r2 + k2, ej_r1 + k1) +
                  BilinearPairingsXInv.blocks[*b].elt(ej_r1 + k1, ej_s2 + k2) *
                  BilinearPairingsY   .blocks[*b].elt(ej_r2 + k2, ej_s1 + k1))/4;
        }
        SchurComplement.blocks[j].elt(u1, u2) = tmp;
        if (u2 != u1)
          SchurComplement.blocks[j].elt(u2, u1) = tmp;
      }
    }
  }
}

// dualResidues[p] = primalObjective[p] - Tr(A_p Y) - (FreeVarMatrix y)_p,
// for 0 <= p < primalObjective.size()
//
// The pairings Tr(A_p Y) can be written in terms of BilinearPairingsY:
//
//   Tr(A_(j,r,s,k) Y) = \sum_{b \in blocks[j]}
//                       (1/2) (BilinearPairingsY_{ej r + k, ej s + k} +
//                              swap (r <-> s))
// where ej = d_j + 1.
//
// Inputs: sdp, y, BilinearPairingsY
// Output: dualResidues (overwriten)
//
void computeDualResidues(const SDP &sdp,
                         const Vector &y,
                         const BlockDiagonalMatrix &BilinearPairingsY,
                         Vector &dualResidues) {
  #pragma omp parallel for schedule(dynamic)
  for (unsigned int j = 0; j < sdp.dimensions.size(); j++) {
    const int ej = sdp.degrees[j] + 1;

    for (vector<IndexTuple>::const_iterator t = sdp.constraintIndices[j].begin();
         t != sdp.constraintIndices[j].end();
         t++) {
      const int p    = t->p;
      const int ej_r = t->r * ej;
      const int ej_s = t->s * ej;
      const int k    = t->k;

      // dualResidues[p] = -Tr(A_p Y)
      dualResidues[p] = 0;
      for (vector<int>::const_iterator b = sdp.blocks[j].begin();
           b != sdp.blocks[j].end(); b++) {
        dualResidues[p] -= BilinearPairingsY.blocks[*b].elt(ej_r+k, ej_s+k);
        dualResidues[p] -= BilinearPairingsY.blocks[*b].elt(ej_s+k, ej_r+k);
      }
      dualResidues[p] /= 2;

      // dualResidues[p] = -Tr(A_p Y) - (FreeVarMatrix y)_p
      for (int n = 0; n < sdp.FreeVarMatrix.cols; n++)
        dualResidues[p] -= sdp.FreeVarMatrix.elt(p, n)*y[n];

      // dualResidues[p] = primalObjective[p] - Tr(A_p Y) - (FreeVarMatrix y)_p
      dualResidues[p] += sdp.primalObjective[p];
    }
  }
}

// Result = \sum_p a[p] A_p,
//
// where a[p] is a vector of length primalObjective.size() and the
// constraint matrices A_p are given by
//
//   A_(j,r,s,k) = \sum_{b \in blocks[j]}
//                     Block_b(v_{b,k} v_{b,k}^T \otimes E^{rs}),
//
// where v_{b,k} is the k-th column of bilinearBases[b], as described
// in SDP.h.
//
// Inputs: sdp, a
// Output: Result (overwritten)
//
void constraintMatrixWeightedSum(const SDP &sdp,
                                 const Vector a,
                                 BlockDiagonalMatrix &Result)  {
  #pragma omp parallel for schedule(dynamic)
  for (unsigned int j = 0; j < sdp.dimensions.size(); j++) {
    const int ej = sdp.degrees[j] + 1;

    // For each j, t points to the first IndexTuple corresponding to j
    for (vector<IndexTuple>::const_iterator t = sdp.constraintIndices[j].begin();
         t != sdp.constraintIndices[j].end();
         t += ej) {
      const int p = t->p;
      const int r = t->r;
      const int s = t->s;
      assert(t->k == 0);

      for (vector<int>::const_iterator b = sdp.blocks[j].begin();
           b != sdp.blocks[j].end(); b++) {

        // Result.blocks[b]^(r,s) = V diag(a') V^T, where
        // V=sdp.bilinearBases[b], a' denotes the subvector of a
        // corresponding to j, and M^(r,s) denotes the (r,s)-th block
        // of M.
        diagonalCongruence(&a[p], sdp.bilinearBases[*b], r, s, Result.blocks[*b]);

        // Result should be symmetric, so if r != s, we must divide
        // the (r,s)-th block of Result.blocks[b] by 2 and copy its
        // transpose to the (s,r)-th block.
        if (r != s) {
          const int u = sdp.bilinearBases[*b].rows;
          for (int m = r*u; m < (r+1)*u; m++) {
            for (int n = s*u; n < (s+1)*u; n++) {
              Result.blocks[*b].elt(m, n) /= 2;
              Result.blocks[*b].elt(n, m) = Result.blocks[*b].elt(m, n);
            }
          }
        }
      }
    }
  }
}

// Compute the vectors r_x and r_y on the right-hand side of the Schur
// complement equation:
//
// {{S, -B}, {B^T, 0}} . {dx, dy} = {r_x, r_y}
//
// where S = SchurComplement and B = FreeVarMatrix.  Specifically,
//
// r_x[p] = -dualResidues[p] - Tr(A_p Z)              for 0 <= p < P
// r_y[n] = dualObjective[n] - (FreeVarMatrix^T x)_n  for 0 <= n < N
//
// where P = primalObjective.size(), N = dualObjective.size()
//
// Inputs:
// - sdp, an SDP
// - dualResidues, a Vector of length P
// - Z = X^{-1} (PrimalResidues Y - R)
// - x, a vector of length P
// Outputs:
// - r_x, a Vector of length P
// - r_y, a Vector of length N
//
void computeSchurRHS(const SDP &sdp,
                     const Vector &dualResidues,
                     const BlockDiagonalMatrix &Z,
                     const Vector &x,
                     Vector &r_x,
                     Vector &r_y) {
  // r_x[p] = -dualResidues[p]
  for (unsigned int p = 0; p < r_x.size(); p++)
    r_x[p] = -dualResidues[p];

  // r_x[p] = -dualResidues[p] - Tr(A_p Z), where A_p are as described
  // in SDP.h.  The trace can be computed in terms of bilinearBases
  // using bilinearBlockPairing.
  #pragma omp parallel for schedule(dynamic)
  for (unsigned int j = 0; j < sdp.dimensions.size(); j++) {
    for (vector<IndexTuple>::const_iterator t = sdp.constraintIndices[j].begin();
         t != sdp.constraintIndices[j].end();
         t++) {
      for (vector<int>::const_iterator b = sdp.blocks[j].begin();
           b != sdp.blocks[j].end(); b++) {
        const int h = sdp.bilinearBases[*b].rows;
        // Pointer to the k-th column of sdp.bilinearBases[*b]
        const Real *q = &sdp.bilinearBases[*b].elements[(t->k) * h];

        r_x[t->p] -= bilinearBlockPairing(q, h, Z.blocks[*b], t->r, t->s);
      }
    }
  }

  // r_y[n] = dualObjective[n] - (FreeVarMatrix^T x)_n
  #pragma omp parallel for schedule(static)
  for (unsigned int n = 0; n < sdp.dualObjective.size(); n++) {
    r_y[n] = sdp.dualObjective[n];
    for (int p = 0; p < sdp.FreeVarMatrix.rows; p++) {
      r_y[n] -= sdp.FreeVarMatrix.elt(p, n)*x[p];
    }
  }
}

// PrimalResidues = \sum_p A_p x[p] - X
//
// Inputs: sdp, x, X
// Output: PrimalResidues (overwritten)
//
void computePrimalResidues(const SDP &sdp,
                           const Vector x,
                           const BlockDiagonalMatrix &X,
                           BlockDiagonalMatrix &PrimalResidues) {
  constraintMatrixWeightedSum(sdp, x, PrimalResidues);
  PrimalResidues -= X;
}

// Centering parameter \beta_p for the predictor step
Real predictorCenteringParameter(const SDPSolverParameters &parameters,
                                 const bool isPrimalDualFeasible) {
  return isPrimalDualFeasible ? Real(0) : parameters.infeasibleCenteringParameter;
}

// Centering parameter \beta_c for the corrector step
Real correctorCenteringParameter(const SDPSolverParameters &parameters,
                                 const BlockDiagonalMatrix &X,
                                 const BlockDiagonalMatrix &dX,
                                 const BlockDiagonalMatrix &Y,
                                 const BlockDiagonalMatrix &dY,
                                 const Real &mu,
                                 const bool isPrimalDualFeasible) {
  Real r = frobeniusProductOfSums(X, dX, Y, dY) / (mu * X.dim);
  Real beta = r < 1 ? r*r : r;

  if (isPrimalDualFeasible)
    return min(max(parameters.feasibleCenteringParameter, beta), Real(1));
  else
    return max(parameters.infeasibleCenteringParameter, beta);
}

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
// - dM, a BlockDiagonalMatrix with the same structure as M
// Workspace:
// - MInvDM (NB: overwritten when computing minEigenvalue)
// - eigenvalues, a Vector of eigenvalues for each block of M
// - workspace, a vector of Vectors needed by the minEigenvalue function
// Output:
// - min(\gamma \alpha(M, dM), 1) (returned)
//
Real stepLength(BlockDiagonalMatrix &MCholesky,
                BlockDiagonalMatrix &dM,
                BlockDiagonalMatrix &MInvDM,
                vector<Vector> &eigenvalues,
                vector<Vector> &workspace,
                const Real gamma) {
  // MInvDM = L^{-1} dM L^{-T}, where M = L L^T
  MInvDM.copyFrom(dM);
  lowerTriangularInverseCongruence(MInvDM, MCholesky);

  const Real lambda = minEigenvalue(MInvDM, workspace, eigenvalues);
  if (lambda > -gamma)
    return 1;
  else
    return -gamma/lambda;
}

// Compute the quantities needed to solve the Schur complement
// equation
//
// {{S, -B}, {B^T, 0}} . {dx, dy} = {r, s}
//
// (where S = SchurComplement, B = FreeVarMatrix), using the method
// described in the manual:
//
// - Compute S using BilinearPairingsXInv and BilinearPairingsY.
//
// - Stabilize S by adding a low-rank update S' = S + U U^T and
//   compute the Cholesky decomposition S' = L' L'^T.
//
// - Form B' = (B U) and compute
//
//   - SchurOffDiagonal = L'^{-1} B
//   - L'^{-1} U (this is stored implicitly in the stabilizeBlock* variables)
//   - Q = (L'^{-1} B')^T (L'^{-1} B') - {{0, 0}, {0, 1}}
//
// - Compute the LU decomposition of Q.
//
// This data is sufficient to efficiently solve the above equation for
// a given r,s.
//
// Inputs:
// - BilinearPairingsXInv, BilinearPairingsY (these are members of
//   SDPSolver, but we include them as arguments to emphasize that
//   they must be computed first)
// - choleskyStabilizeThreshold: the real constant \theta used in
//   stabilizing S
// Workspace (members of SDPSolver which are modified by this method
// and not used later):
// - SchurComplement
// - schurStabilizeIndices
// - schurStabilizeLambdas
// Outputs (members of SDPSolver which are modified by this method and
// used later):
// - SchurComplementCholesky
// - SchurOffDiagonal
// - stabilizeBlockIndices
// - stabilizeBlockUpdateRow
// - stabilizeBlockUpdateColumn
// - stabilizeBlocks
// - Q, Qpivots
//
void SDPSolver::initializeSchurComplementSolver(const BlockDiagonalMatrix &BilinearPairingsXInv,
                                                const BlockDiagonalMatrix &BilinearPairingsY,
                                                const Real &choleskyStabilizeThreshold) {
  computeSchurComplement(sdp, BilinearPairingsXInv, BilinearPairingsY, SchurComplement);

  // compute SchurComplementCholesky = L', where
  //
  //   L' L'^T = S' = SchurComplement + U U^T
  //
  // Here, the 'update' matrix U has columns given by
  //
  //   U = ( Lambda_{p_1} e_{p_1}, ..., Lambda_{p_M} e_{p_M} )
  //
  // where e_p is a unit vector in the p-th direction and the
  // Lambda_{p_m} are constants. The p_i are given by
  // schurStabilizeIndices and the corresponding Lambda_i by
  // schurStabilizeLambdas.
  // 
  choleskyDecompositionStabilized(SchurComplement, SchurComplementCholesky,
                                  schurStabilizeIndices,
                                  schurStabilizeLambdas,
                                  cast2double(choleskyStabilizeThreshold));

  // SchurOffDiagonal = L'^{-1} FreeVarMatrix
  SchurOffDiagonal.copyFrom(sdp.FreeVarMatrix);
  blockMatrixLowerTriangularSolve(SchurComplementCholesky, SchurOffDiagonal);

  // Next we compute L'^{-1} U, which is stored implicitly in terms of
  // its nonzero submatrices in the stabilizeBlock* variables.

  // total number of columns in the off-diagonal part B' = (B U)
  // (currently just B; will accumulate the rest shortly)
  int offDiagonalColumns = SchurOffDiagonal.cols;

  stabilizeBlockIndices.resize(0);
  stabilizeBlockUpdateRow.resize(0);
  stabilizeBlockUpdateColumn.resize(0);

  // j runs over blocks of SchurComplement
  for (unsigned int j = 0; j < SchurComplement.blocks.size(); j++) {
    if (schurStabilizeIndices[j].size() > 0) {
      // the j-th block of S contains stabilized directions. We have a
      // block submatrix
      //
      //   U_j = stabilizeBlocks[j]
      //
      // of U for each such j.
      stabilizeBlockIndices.push_back(j);
      
      // index of the first stabilized direction within the j-th block
      int startIndex = schurStabilizeIndices[j][0];

      // set the dimensions of U_j 
      stabilizeBlocks[j].resize(SchurComplement.blocks[j].rows - startIndex,
                                schurStabilizeIndices[j].size());
      // set U_j = 0
      stabilizeBlocks[j].setZero();
      // for each column of U_j add Lambda_p in the row (p - startIndex)
      for (unsigned int c = 0; c < schurStabilizeIndices[j].size(); c++) {
        int r = schurStabilizeIndices[j][c] - startIndex;
        stabilizeBlocks[j].elt(r, c) = schurStabilizeLambdas[j][c];
      }

      // append the row of U corresponding to the top-left of U_j
      stabilizeBlockUpdateRow.push_back(SchurComplement.blockStartIndices[j] + startIndex);
      // append the column of U corresponding to the top-left of U_j
      stabilizeBlockUpdateColumn.push_back(offDiagonalColumns);
      // update the number of off-diagonal columns
      offDiagonalColumns += stabilizeBlocks[j].cols;
    }
  }

  // Set U = L'^{-1} U
  //
  // We do this by modifying the blocks U_j = stabilizeBlocks[j]
  // in-place, multiplying by the inverse of the appropriate submatrix
  // of SchurComplementCholesky.  We henceforth refer to L'^{-1} U as
  // V to avoid confusion.
  #pragma omp parallel for schedule(dynamic)
  for (unsigned int j = 0; j < stabilizeBlockIndices.size(); j++) {
    int b = stabilizeBlockIndices[j];
    int startIndex = schurStabilizeIndices[b][0];
    Rtrsm("Left", "Lower", "NoTranspose", "NonUnitDiagonal",
          stabilizeBlocks[b].rows, stabilizeBlocks[b].cols, 1,
          &SchurComplementCholesky.blocks[b].elt(startIndex, startIndex),
          SchurComplementCholesky.blocks[b].rows,
          &stabilizeBlocks[b].elt(0, 0),
          stabilizeBlocks[b].rows);
  }

  // Next, we compute
  //
  //   Q = (L'^{-1} B')^T (L'^{-1} B') - {{0, 0}, {0, 1}}
  //
  // Where B' = (B U).  We think of Q as containing four blocks called
  // Upper/Lower-Left/Right.

  // Set the dimensions of Q
  Q.resize(offDiagonalColumns, offDiagonalColumns);
  Q.setZero();

  // Here, SchurOffDiagonal = L'^{-1} B.
  //
  // UpperLeft(Q) = SchurOffDiagonal^T SchurOffDiagonal
  matrixSquareIntoBlock(SchurOffDiagonal, Q, 0, 0);

  // Here, stabilizeBlocks contains the blocks of V = L'^{-1} U.
  //
  // LowerRight(Q) = V^T V - 1
  for (unsigned int j = 0; j < stabilizeBlockIndices.size(); j++) {
    int b = stabilizeBlockIndices[j];
    int c = stabilizeBlockUpdateColumn[j];
    matrixSquareIntoBlock(stabilizeBlocks[b], Q, c, c);
    // subtract the identity matrix from this block
    for (int i = c; i < c + stabilizeBlocks[b].cols; i++)
      Q.elt(i, i) -= 1;
  }

  // LowerLeft(Q) = V^T SchurOffDiagonal
  # pragma omp parallel for schedule(dynamic)
  for (unsigned int j = 0; j < stabilizeBlockIndices.size(); j++) {
    int b = stabilizeBlockIndices[j];
    int p = stabilizeBlockUpdateRow[j];
    int r = stabilizeBlockUpdateColumn[j];
    Rgemm("Transpose", "NoTranspose",
          stabilizeBlocks[b].cols,
          SchurOffDiagonal.cols,
          stabilizeBlocks[b].rows,
          1,
          &stabilizeBlocks[b].elements[0],
          stabilizeBlocks[b].rows,
          &SchurOffDiagonal.elt(p, 0),
          SchurOffDiagonal.rows,
          0,
          &Q.elt(r, 0),
          Q.rows);
  }

  // UpperRight(Q) = LowerLeft(Q)^T
  # pragma omp parallel for schedule(static)
  for (int c = 0; c < SchurOffDiagonal.cols; c++)
    for (int r = SchurOffDiagonal.cols; r < Q.rows; r++)
      Q.elt(c, r) = Q.elt(r, c);

  // Resize Qpivots appropriately and compute the LU decomposition of Q
  Qpivots.resize(Q.rows);
  LUDecomposition(Q, Qpivots);
}


// Solve the Schur complement equation for dx, dy.
//
// - As inputs, dx and dy are the residues r_x and r_y on the
//   right-hand side of the Schur complement equation.
// - As outputs, dx and dy are overwritten with the solutions of the
//   Schur complement equation.
//
// The equation is solved using the block-decomposition described in
// the manual.
//
void SDPSolver::solveSchurComplementEquation(Vector &dx, Vector &dy) {
  // dx = SchurComplementCholesky^{-1} dx
  blockMatrixLowerTriangularSolve(SchurComplementCholesky, dx);

  // extend dy by additional coordinates z needed for stabilizing
  // dyExtended = (dy, z)
  dyExtended.resize(Q.rows);

  // dy = r_y - SchurOffDiagonal^T dx
  for (unsigned int n = 0; n < dy.size(); n++)
    dyExtended[n] = dy[n];
  vectorScaleMatrixMultiplyTransposeAdd(-1, SchurOffDiagonal, dx, 1, dyExtended);

  // z = -V^T dx
  for (unsigned int j = 0; j < stabilizeBlockIndices.size(); j++) {
    int b = stabilizeBlockIndices[j];
    int pTopLeft = stabilizeBlockUpdateRow[j];
    int cTopLeft = stabilizeBlockUpdateColumn[j];

    for (int c = 0; c < stabilizeBlocks[b].cols; c++) {
      dyExtended[cTopLeft + c] = 0;
      for (int r = 0; r < stabilizeBlocks[b].rows; r++)
        dyExtended[cTopLeft + c] -= dx[pTopLeft + r] * stabilizeBlocks[b].elt(r, c);
    }
  }

  // dyExtended = Q^{-1} dyExtended
  solveWithLUDecomposition(Q, Qpivots, dyExtended);

  // dx += SchurOffDiagonal dy
  vectorScaleMatrixMultiplyAdd(1, SchurOffDiagonal, dyExtended, 1, dx);

  // dx += V z
  for (unsigned int j = 0; j < stabilizeBlockIndices.size(); j++) {
    int b = stabilizeBlockIndices[j];
    int pTopLeft = stabilizeBlockUpdateRow[j];
    int cTopLeft = stabilizeBlockUpdateColumn[j];

    for (int c = 0; c < stabilizeBlocks[b].cols; c++)
      for (int r = 0; r < stabilizeBlocks[b].rows; r++)
        dx[pTopLeft + r] += dyExtended[cTopLeft + c] * stabilizeBlocks[b].elt(r, c);
  }

  // dx = SchurComplementCholesky^{-T} dx
  blockMatrixLowerTriangularTransposeSolve(SchurComplementCholesky, dx);

  // dy = first few entries of dyExtended
  for (unsigned int n = 0; n < dy.size(); n++)
    dy[n] = dyExtended[n];
}

// Compute the search direction (dx, dX, dy, dY) for the predictor and
// corrector phases.
//
// Inputs:
// - beta, the centering parameter
// - mu = Tr(X Y) / X.cols
// - correctorPhase: boolean indicating whether we're in the corrector
//   phase or predictor phase.
// Workspace (members of SDPSolver which are modified in-place but not
// used elsewhere):
// - Z, R
// Outputs (members of SDPSolver which are modified in-place):
// - dx, dX, dy, dY
//
void SDPSolver::computeSearchDirection(const Real &beta,
                                       const Real &mu,
                                       const bool correctorPhase) {
  // R = beta mu I - X Y (predictor phase)
  // R = beta mu I - X Y - dX dY (corrector phase)
  blockDiagonalMatrixScaleMultiplyAdd(-1, X, Y, 0, R);
  if (correctorPhase)
    blockDiagonalMatrixScaleMultiplyAdd(-1, dX, dY, 1, R);
  R.addDiagonal(beta*mu);

  // Z = Symmetrize(X^{-1} (PrimalResidues Y - R))
  blockDiagonalMatrixMultiply(PrimalResidues, Y, Z);
  Z -= R;
  blockMatrixSolveWithCholesky(XCholesky, Z);
  Z.symmetrize();

  // r_x[p] = -dualResidues[p] - Tr(A_p Z)
  // r_y[n] = dualObjective[n] - (FreeVarMatrix^T x)_n
  // Here, dx = r_x, dy = r_y.
  computeSchurRHS(sdp, dualResidues, Z, x, dx, dy);

  // Solve for dx, dy in-place
  solveSchurComplementEquation(dx, dy);

  // dX = PrimalResidues + \sum_p A_p dx[p]
  constraintMatrixWeightedSum(sdp, dx, dX);
  dX += PrimalResidues;

  // dY = Symmetrize(X^{-1} (R - dX Y))
  blockDiagonalMatrixMultiply(dX, Y, dY);
  dY -= R;
  blockMatrixSolveWithCholesky(XCholesky, dY);
  dY.symmetrize();
  dY *= -1;
}

/***********************************************************************/
// The main solver loop

SDPSolverTerminateReason SDPSolver::run(const path checkpointFile) {
  Real primalStepLength;
  Real dualStepLength;

  printHeader();

  for (int iteration = 1;; iteration++) {
    if (timers["Last checkpoint"].elapsed().wall >= parameters.checkpointInterval * 1000000000LL)
      saveCheckpoint(checkpointFile);
    if (timers["Solver runtime"].elapsed().wall >= parameters.maxRuntime * 1000000000LL)
      return MaxRuntimeExceeded;

    primalObjective = sdp.objectiveConst + dotProduct(sdp.primalObjective, x);
    dualObjective   = sdp.objectiveConst + dotProduct(sdp.dualObjective, y);
    dualityGap      = abs(primalObjective - dualObjective) /
      max(Real(abs(primalObjective) + abs(dualObjective)), Real(1));

    choleskyDecomposition(X, XCholesky);
    choleskyDecomposition(Y, YCholesky);

    // Compute the bilinear pairings BilinearPairingsXInv and
    // BilinearPairingsY needed for the dualResidues and the Schur
    // complement matrix
    blockTensorInvTransposeCongruenceWithCholesky(XCholesky, sdp.bilinearBases,
                                                  bilinearPairingsWorkspace,
                                                  BilinearPairingsXInv);
    blockTensorTransposeCongruence(Y, sdp.bilinearBases,
                                   bilinearPairingsWorkspace,
                                   BilinearPairingsY);

    // dualResidues[p] = primalObjective[p] - Tr(A_p Y) - (FreeVarMatrix y)_p,
    computeDualResidues(sdp, y, BilinearPairingsY, dualResidues);
    dualError = maxAbsVector(dualResidues);

    // PrimalResidues = \sum_p A_p x[p] - X
    computePrimalResidues(sdp, x, X, PrimalResidues);
    primalError = PrimalResidues.maxAbs();

    const bool isPrimalFeasible = primalError < parameters.primalErrorThreshold;
    const bool isDualFeasible   = dualError   < parameters.dualErrorThreshold;
    const bool isOptimal        = dualityGap  < parameters.dualityGapThreshold;

    if (isPrimalFeasible && isDualFeasible && isOptimal)                   return PrimalDualOptimal;
    else if (isPrimalFeasible && parameters.findPrimalFeasible)            return PrimalFeasible;
    else if (isDualFeasible && parameters.findDualFeasible)                return DualFeasible;
    else if (primalStepLength == 1 && parameters.detectPrimalFeasibleJump) return PrimalFeasibleJumpDetected;
    else if (dualStepLength == 1 && parameters.detectDualFeasibleJump)     return DualFeasibleJumpDetected;
    else if (iteration > parameters.maxIterations)                         return MaxIterationsExceeded;

    // Compute SchurComplement and prepare to solve the Schur
    // complement equation for dx, dy
    initializeSchurComplementSolver(BilinearPairingsXInv, BilinearPairingsY,
                                    parameters.choleskyStabilizeThreshold);

    // Compute the complementarity mu = Tr(X Y)/X.dim
    Real mu = frobeniusProductSymmetric(X, Y)/X.dim;
    if (mu > parameters.maxComplementarity)
      return MaxComplementarityExceeded;

    // Compute the predictor solution for (dx, dX, dy, dY)
    Real betaPredictor =
      predictorCenteringParameter(parameters, isPrimalFeasible && isDualFeasible);
    computeSearchDirection(betaPredictor, mu, false);

    // Compute the corrector solution for (dx, dX, dy, dY)
    Real betaCorrector =
      correctorCenteringParameter(parameters, X, dX, Y, dY, mu,
                                  isPrimalFeasible && isDualFeasible);
    computeSearchDirection(betaCorrector, mu, true);

    // Compute step-lengths that preserve positive definiteness of X, Y
    primalStepLength = stepLength(XCholesky, dX, StepMatrixWorkspace,
                                  eigenvaluesWorkspace, QRWorkspace,
                                  parameters.stepLengthReduction);
    dualStepLength   = stepLength(YCholesky, dY, StepMatrixWorkspace,
                                  eigenvaluesWorkspace, QRWorkspace,
                                  parameters.stepLengthReduction);

    // If our problem is both dual-feasible and primal-feasible,
    // ensure we're following the true Newton direction.
    if (isPrimalFeasible && isDualFeasible) {
      primalStepLength = min(primalStepLength, dualStepLength);
      dualStepLength = primalStepLength;
    }

    printIteration(iteration, mu, primalStepLength, dualStepLength, betaCorrector);

    // Update the primal point (x, X) += primalStepLength*(dx, dX)
    addScaledVector(x, primalStepLength, dx);
    dX *= primalStepLength;
    X += dX;

    // Update the dual point (y, Y) += dualStepLength*(dy, dY)
    addScaledVector(y, dualStepLength, dy);
    dY *= dualStepLength;
    Y += dY;
  }

  // Never reached
  return MaxIterationsExceeded;
}
