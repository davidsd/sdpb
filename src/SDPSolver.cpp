//=======================================================================
// Copyright 2014 David Simmons-Duffin.
// Distributed under the MIT License.
// (See accompanying file LICENSE or copy at
//  http://opensource.org/licenses/MIT)
//=======================================================================


#include <algorithm>
#include <iostream>
#include <vector>
#include "omp.h"
#include "boost/filesystem.hpp"
#include "SDPSolver.h"
#include "Timers.h"

using boost::filesystem::path;
using boost::timer::nanosecond_type;
using std::cout;

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
  SchurBlocks(sdp.schurBlockDims()),
  SchurBlocksCholesky(SchurBlocks),
  SchurOffDiagonal(sdp.FreeVarMatrix),
  schurStabilizeIndices(SchurBlocks.blocks.size()),
  schurStabilizeLambdas(SchurBlocks.blocks.size()),
  stabilizeBlocks(SchurBlocks.blocks.size()),
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
// - Work   : l*m x n*m matrix, intermediate workspace (overwritten)
// - Result : n*m x n*m symmetric matrix (overwritten)
//
// An explanation of the name: a 'congruence' refers to the action M
// -> A M A^T.  We use 'congruence transpose' to refer to a congruence
// with the transposed matrix M -> A^T M A.  'tensor' here refers to
// the fact that we're performing a congruence with the tensor product
// Q \otimes 1.
//
void tensorMatrixCongruenceTranspose(const Matrix &A,
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

// Result = Q'^T A^{-1} Q', where Q' = Q \otimes 1, where \otimes
// denotes tensor product and 1 is an mxm idenity matrix.
// Inputs:
// - L      : l*m x l*m cholesky decomposition of A
// - Q      : l   x n   matrix
// - Work   : l*m x n*m matrix, intermediate workspace (overwritten)
// - Result : n*m x n*m symmetric matrix (overwritten)
//
void tensorMatrixInvCongruenceTransposeWithCholesky(const Matrix &L,
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

// v^T A' v, where A' is the (blockRow,blockCol)-th dim x dim block
// inside A.
//
// Input:
// - v        : pointer to the beginning of a vector of length dim
// - dim      : length of the vector v
// - A        : (k*dim) x (k*dim) matrix, where k > blockRow, blockCol
// - blockRow : integer labeling block of A
// - blockCol : integer labeling block of A
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
// Subroutines needed for each iteration

// Result_b = Q[b]'^T A_b Q[b]' for each block 0 <= b < Q.size()
// - Result_b, A_b denote the b-th blocks of Result, A, resp.
// - Q[b]' = Q[b] \otimes 1, where \otimes denotes tensor product
//
// This is just tensorMatrixCongruenceTranspose for each block of a
// BlockDiagonalMatrix.
//
void computeBilinearPairings(const BlockDiagonalMatrix &A,
                             const vector<Matrix> &Q,
                             vector<Matrix> &workspace,
                             BlockDiagonalMatrix &Result) {
  #pragma omp parallel for schedule(dynamic)
  for (unsigned int b = 0; b < Q.size(); b++)
    tensorMatrixCongruenceTranspose(A.blocks[b], Q[b],
                                    workspace[b], Result.blocks[b]);
}

// Result_b = Q[b]'^T A_b^{-1} Q[b]' for each block 0 <= b < Q.size()
// - Result_b, A_b denote the b-th blocks of Result, A, resp.
// - Q[b]' = Q[b] \otimes 1, where \otimes denotes tensor product
// 
// This is just tensorMatrixInvCongruenceTransposeWithCholesky for
// each block of a BlockDiagonalMatrix.
//
void computeInvBilinearPairingsWithCholesky(const BlockDiagonalMatrix &L,
                                            const vector<Matrix> &Q,
                                            vector<Matrix> &workspace,
                                            BlockDiagonalMatrix &Result) {
  #pragma omp parallel for schedule(dynamic)
  for (unsigned int b = 0; b < Q.size(); b++)
    tensorMatrixInvCongruenceTransposeWithCholesky(L.blocks[b],
                                                   Q[b],
                                                   workspace[b],
                                                   Result.blocks[b]);
}

void computeSchurBlocks(const SDP &sdp,
                        const BlockDiagonalMatrix &BilinearPairingsXInv,
                        const BlockDiagonalMatrix &BilinearPairingsY,
                        BlockDiagonalMatrix &SchurBlocks) {
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
        SchurBlocks.blocks[j].elt(u1, u2) = tmp;
        if (u2 != u1)
          SchurBlocks.blocks[j].elt(u2, u1) = tmp;
      }
    }
  }
}

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

      dualResidues[p] = 0;
      for (vector<int>::const_iterator b = sdp.blocks[j].begin();
           b != sdp.blocks[j].end(); b++) {
        dualResidues[p] -= BilinearPairingsY.blocks[*b].elt(ej_r+k, ej_s+k);
        dualResidues[p] -= BilinearPairingsY.blocks[*b].elt(ej_s+k, ej_r+k);
      }
      dualResidues[p] /= 2;

      for (int n = 0; n < sdp.FreeVarMatrix.cols; n++)
        dualResidues[p] -= sdp.FreeVarMatrix.elt(p, n)*y[n];
      dualResidues[p] += sdp.primalObjective[p];
    }
  }
}

void constraintMatrixWeightedSum(const SDP &sdp,
                                 const Vector x,
                                 BlockDiagonalMatrix &result)  {
  #pragma omp parallel for schedule(dynamic)
  for (unsigned int j = 0; j < sdp.dimensions.size(); j++) {
    const int ej = sdp.degrees[j] + 1;

    for (vector<IndexTuple>::const_iterator t = sdp.constraintIndices[j].begin();
         t != sdp.constraintIndices[j].end();
         t += ej) {
      const int p = t->p;
      const int r = t->r;
      const int s = t->s;
      assert(t->k == 0);

      for (vector<int>::const_iterator b = sdp.blocks[j].begin();
           b != sdp.blocks[j].end(); b++) {
        diagonalCongruence(&x[p], sdp.bilinearBases[*b], r, s,
                           result.blocks[*b]);
        if (r != s) {
          const int u = sdp.bilinearBases[*b].rows;
          for (int m = r*u; m < (r+1)*u; m++) {
            for (int n = s*u; n < (s+1)*u; n++) {
              result.blocks[*b].elt(m, n) /= 2;
              result.blocks[*b].elt(n, m) = result.blocks[*b].elt(m, n);
            }
          }
        }
      }
    }
  }
}

void computeSchurRHS(const SDP &sdp,
                     const Vector &dualResidues,
                     const BlockDiagonalMatrix &Z,
                     const Vector &x,
                     Vector &r,
                     Vector &s) {
  for (unsigned int p = 0; p < r.size(); p++)
    r[p] = -dualResidues[p];

  #pragma omp parallel for schedule(dynamic)
  for (unsigned int j = 0; j < sdp.dimensions.size(); j++) {
    for (vector<IndexTuple>::const_iterator t = sdp.constraintIndices[j].begin();
         t != sdp.constraintIndices[j].end();
         t++) {
      for (vector<int>::const_iterator b = sdp.blocks[j].begin();
           b != sdp.blocks[j].end(); b++) {
        const int delta = sdp.bilinearBases[*b].rows;
        // Pointer to the k-th column of sdp.bilinearBases[*b]
        const Real *q = &sdp.bilinearBases[*b].elements[(t->k) * delta];

        r[t->p] -= bilinearBlockPairing(q, delta, Z.blocks[*b], t->r, t->s);
      }
    }
  }

  #pragma omp parallel for schedule(static)
  for (unsigned int n = 0; n < sdp.dualObjective.size(); n++) {
    s[n] = sdp.dualObjective[n];
    for (int p = 0; p < sdp.FreeVarMatrix.rows; p++) {
      s[n] -= sdp.FreeVarMatrix.elt(p, n)*x[p];
    }
  }
}

// PrimalResidues = sum_p F_p x_p - X - F_0
//
void computePrimalResidues(const SDP &sdp,
                           const Vector x,
                           const BlockDiagonalMatrix &X,
                           BlockDiagonalMatrix &PrimalResidues) {
  constraintMatrixWeightedSum(sdp, x, PrimalResidues);
  PrimalResidues -= X;
}

Real predictorCenteringParameter(const SDPSolverParameters &parameters,
                                 const bool isPrimalDualFeasible) {
  return isPrimalDualFeasible ? Real(0) : parameters.infeasibleCenteringParameter;
}

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

Real stepLength(BlockDiagonalMatrix &XCholesky,
                BlockDiagonalMatrix &dX,
                BlockDiagonalMatrix &XInvDX,
                vector<Vector> &eigenvalues,
                vector<Vector> &workspace,
                const SDPSolverParameters &parameters) {
  // XInvDX = L^{-1} dX L^{-1 T}, where X = L L^T
  XInvDX.copyFrom(dX);
  lowerTriangularInverseCongruence(XInvDX, XCholesky);

  const Real lambda = minEigenvalue(XInvDX, workspace, eigenvalues);
  const Real gamma  = parameters.stepLengthReduction;
  if (lambda > -gamma)
    return 1;
  else
    return -gamma/lambda;
}

void SDPSolver::initializeSchurComplementSolver(const BlockDiagonalMatrix &BilinearPairingsXInv,
                                                const BlockDiagonalMatrix &BilinearPairingsY,
                                                const Real &choleskyStabilizeThreshold) {
  computeSchurBlocks(sdp, BilinearPairingsXInv, BilinearPairingsY, SchurBlocks);
  choleskyDecompositionStabilized(SchurBlocks, SchurBlocksCholesky,
                                  schurStabilizeIndices,
                                  schurStabilizeLambdas,
                                  cast2double(choleskyStabilizeThreshold));

  // SchurOffDiagonal = {{- 1, 0}, {E, G}}
  SchurOffDiagonal.copyFrom(sdp.FreeVarMatrix);
  blockMatrixLowerTriangularSolve(SchurBlocksCholesky, SchurOffDiagonal);
  int updateColumns = SchurOffDiagonal.cols;

  stabilizeBlockIndices.resize(0);
  stabilizeBlockUpdateRow.resize(0);
  stabilizeBlockUpdateColumn.resize(0);

  for (unsigned int j = 0; j < SchurBlocks.blocks.size(); j++) {
    if (schurStabilizeIndices[j].size() > 0) {
      int startIndex = schurStabilizeIndices[j][0];
      int blockRows  = SchurBlocks.blocks[j].rows - startIndex;
      int blockCols  = schurStabilizeIndices[j].size();

      stabilizeBlockIndices.push_back(j);
      stabilizeBlockUpdateRow.push_back(SchurBlocks.blockStartIndices[j] + startIndex);
      stabilizeBlockUpdateColumn.push_back(updateColumns);
      updateColumns += blockCols;

      stabilizeBlocks[j].setRowsCols(blockRows, blockCols);
      stabilizeBlocks[j].setZero();
      for (unsigned int c = 0; c < schurStabilizeIndices[j].size(); c++) {
        int r = schurStabilizeIndices[j][c] - startIndex;
        stabilizeBlocks[j].elt(r, c) = schurStabilizeLambdas[j][c];
      }
    }
  }

  // G := SchurBlocksCholesky^{-1} G
  #pragma omp parallel for schedule(dynamic)
  for (unsigned int j = 0; j < stabilizeBlockIndices.size(); j++) {
    int b = stabilizeBlockIndices[j];
    int startIndex = schurStabilizeIndices[b][0];
    Rtrsm("Left", "Lower", "NoTranspose", "NonUnitDiagonal",
          stabilizeBlocks[b].rows, stabilizeBlocks[b].cols, 1,
          &SchurBlocksCholesky.blocks[b].elt(startIndex, startIndex),
          SchurBlocksCholesky.blocks[b].rows,
          &stabilizeBlocks[b].elt(0, 0),
          stabilizeBlocks[b].rows);
  }

  // Q = SchurOffDiagonal^T SchurOffDiagonal - {{0,0},{0,1}}
  Q.setRowsCols(updateColumns, updateColumns);
  Q.setZero();

  matrixSquareIntoBlock(SchurOffDiagonal, Q, 0, 0);

  // LowerRight(Q) = G^T G - 1
  for (unsigned int j = 0; j < stabilizeBlockIndices.size(); j++) {
    int b = stabilizeBlockIndices[j];
    int c = stabilizeBlockUpdateColumn[j];
    matrixSquareIntoBlock(stabilizeBlocks[b], Q, c, c);
    for (int i = c; i < c + stabilizeBlocks[b].cols; i++)
      Q.elt(i, i) -= 1;
  }

  // LowerLeft(Q) = G^T U
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

  Qpivots.resize(Q.rows);
  LUDecomposition(Q, Qpivots);
}


// As inputs, dx and dy are the residues r_x and r_y on the right-hand
// side of the Schur complement equation. As outputs, they are the
// values for dx and dy.
//
void SDPSolver::solveSchurComplementEquation(Vector &dx, Vector &dy) {
  // dx = SchurBlocksCholesky^{-1} dx
  blockMatrixLowerTriangularSolve(SchurBlocksCholesky, dx);

  dyExtended.resize(Q.rows);
  // k_1 = r_y - SchurOffDiagonal^T dx
  for (unsigned int n = 0; n < dy.size(); n++)
    dyExtended[n] = dy[n];

  vectorScaleMatrixMultiplyTransposeAdd(-1, SchurOffDiagonal, dx, 1, dyExtended);

  // k_2 = -G^T dx
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

  // k = Q^{-1} k
  solveWithLUDecomposition(Q, Qpivots, dyExtended);

  // dx = dx + SchurOffDiagonal k_1
  vectorScaleMatrixMultiplyAdd(1, SchurOffDiagonal, dyExtended, 1, dx);
  // dx += G k_2
  for (unsigned int j = 0; j < stabilizeBlockIndices.size(); j++) {
    int b = stabilizeBlockIndices[j];
    int pTopLeft = stabilizeBlockUpdateRow[j];
    int cTopLeft = stabilizeBlockUpdateColumn[j];

    for (int c = 0; c < stabilizeBlocks[b].cols; c++)
      for (int r = 0; r < stabilizeBlocks[b].rows; r++)
        dx[pTopLeft + r] += dyExtended[cTopLeft + c] * stabilizeBlocks[b].elt(r, c);
  }

  // dx = SchurBlocksCholesky^{-T} dx
  blockMatrixLowerTriangularTransposeSolve(SchurBlocksCholesky, dx);
  // dy = k_1
  for (unsigned int n = 0; n < dy.size(); n++)
    dy[n] = dyExtended[n];
}

void SDPSolver::computeSearchDirection(const Real &beta,
                                       const Real &mu,
                                       const bool correctorPhase) {
  blockDiagonalMatrixScaleMultiplyAdd(-1, X, Y, 0, R);
  if (correctorPhase)
    blockDiagonalMatrixScaleMultiplyAdd(-1, dX, dY, 1, R);
  R.addDiagonal(beta*mu);

  // Z = Symmetrize(X^{-1} (PrimalResidues Y - R))
  blockDiagonalMatrixMultiply(PrimalResidues, Y, Z);
  Z -= R;
  blockMatrixSolveWithCholesky(XCholesky, Z);
  Z.symmetrize();

  // (dx_RHS)_k = -d_k + Tr(F_k Z)
  // dz_RHS = f - D^T x
  computeSchurRHS(sdp, dualResidues, Z, x, dx, dy);

  // dx_N = B_{NN}^{-1} dx_N, dx_B = E^T dx_N
  solveSchurComplementEquation(dx, dy);

  // dX = R_p + sum_p F_p dx_p
  constraintMatrixWeightedSum(sdp, dx, dX);
  dX += PrimalResidues;

  // dY = Symmetrize(X^{-1} (R - dX Y))
  blockDiagonalMatrixMultiply(dX, Y, dY);
  dY -= R;
  blockMatrixSolveWithCholesky(XCholesky, dY);
  dY.symmetrize();
  dY *= -1;
}

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

    computeInvBilinearPairingsWithCholesky(XCholesky, sdp.bilinearBases,
                                           bilinearPairingsWorkspace,
                                           BilinearPairingsXInv);
    computeBilinearPairings(Y, sdp.bilinearBases,
                            bilinearPairingsWorkspace,
                            BilinearPairingsY);

    // d_k = c_k - Tr(F_k Y) - (D y)_k
    computeDualResidues(sdp, y, BilinearPairingsY, dualResidues);
    dualError = maxAbsVector(dualResidues);

    // PrimalResidues = sum_p F_p x_p - X - F_0 (F_0 is zero for now)
    computePrimalResidues(sdp, x, X, PrimalResidues);
    primalError = PrimalResidues.maxAbs();

    const bool isPrimalFeasible = primalError < parameters.primalErrorThreshold;
    const bool isDualFeasible   = dualError   < parameters.dualErrorThreshold;
    const bool isOptimal        = dualityGap  < parameters.dualityGapThreshold;

    if (isPrimalFeasible && isDualFeasible && isOptimal)
      return PrimalDualOptimal;
    else if (isPrimalFeasible && parameters.findPrimalFeasible)
      return PrimalFeasible;
    else if (isDualFeasible && parameters.findDualFeasible)
      return DualFeasible;
    else if (primalStepLength == 1 && parameters.detectPrimalFeasibleJump)
      return PrimalFeasibleJumpDetected;
    else if (dualStepLength == 1 && parameters.detectDualFeasibleJump)
      return DualFeasibleJumpDetected;
    // Detect max iterations after computing errors and objective
    // functions for the current point
    else if (iteration > parameters.maxIterations)
      return MaxIterationsExceeded;

    initializeSchurComplementSolver(BilinearPairingsXInv, BilinearPairingsY,
                                    parameters.choleskyStabilizeThreshold);

    Real mu = frobeniusProductSymmetric(X, Y)/X.dim;
    if (mu > parameters.maxComplementarity)
      return MaxComplementarityExceeded;

    // Mehrotra predictor solution for (dx, dX, dY)
    Real betaPredictor = predictorCenteringParameter(parameters,
                                                     isPrimalFeasible && isDualFeasible);
    computeSearchDirection(betaPredictor, mu, false);

    // Mehrotra corrector solution for (dx, dX, dY)
    Real betaCorrector = correctorCenteringParameter(parameters, X, dX, Y, dY, mu,
                                                     isPrimalFeasible && isDualFeasible);
    computeSearchDirection(betaCorrector, mu, true);

    // Step length to preserve positive definiteness
    primalStepLength = stepLength(XCholesky, dX, StepMatrixWorkspace,
                                  eigenvaluesWorkspace, QRWorkspace, parameters);
    dualStepLength   = stepLength(YCholesky, dY, StepMatrixWorkspace,
                                  eigenvaluesWorkspace, QRWorkspace, parameters);

    if (isPrimalFeasible && isDualFeasible) {
      primalStepLength = min(primalStepLength, dualStepLength);
      dualStepLength = primalStepLength;
    }

    printIteration(iteration, mu, primalStepLength, dualStepLength, betaCorrector);

    // Update current point
    scaleVector(dx, primalStepLength);
    addVector(x, dx);
    dX *= primalStepLength;
    X += dX;
    scaleVector(dy, dualStepLength);
    addVector(y, dy);
    dY *= dualStepLength;
    Y += dY;
  }

  // Never reached
  return MaxIterationsExceeded;
}
