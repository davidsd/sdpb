#include <iostream>
#include "omp.h"
#include "boost/filesystem.hpp"
#include "SDPSolver.h"
#include "Timers.h"

using boost::filesystem::path;
using boost::timer::nanosecond_type;
using std::cout;

SDPSolver::SDPSolver(const SDP &sdp):
  sdp(sdp),
  x(sdp.primalObjective.size(), 0),
  X(sdp.psdMatrixBlockDims()),
  y(sdp.dualObjective.size(), 0),
  Y(X),
  dx(x),
  dX(X),
  dy(y),
  dY(Y),
  dualResidues(x),
  PrimalResidues(X),
  XCholesky(X),
  YCholesky(X),
  Z(X),
  R(X),
  BilinearPairingsXInv(sdp.bilinearPairingBlockDims()),
  BilinearPairingsY(BilinearPairingsXInv),
  SchurBlocks(sdp.schurBlockDims()),
  SchurBlocksCholesky(SchurBlocks),
  SchurUpdateLowRank(sdp.FreeVarMatrix),
  schurStabilizeIndices(SchurBlocks.blocks.size()),
  schurStabilizeLambdas(SchurBlocks.blocks.size()),
  stabilizeBlocks(SchurBlocks.blocks.size()),
  Q(sdp.FreeVarMatrix.cols, sdp.FreeVarMatrix.cols),
  Qpivots(sdp.FreeVarMatrix.cols),
  basicKernelCoords(Q.rows),
  StepMatrixWorkspace(X)
{
  // initialize bilinearPairingsWorkspace, eigenvaluesWorkspace, QRWorkspace 
  for (unsigned int b = 0; b < sdp.bilinearBases.size(); b++) {
    bilinearPairingsWorkspace.push_back(Matrix(X.blocks[b].rows, BilinearPairingsXInv.blocks[b].cols));
    eigenvaluesWorkspace.push_back(Vector(X.blocks[b].rows));
    QRWorkspace.push_back(Vector(3*X.blocks[b].rows - 1));
  }
}

// result = b'^T a b', where b' = b \otimes 1
// Inputs:
// - a      : l*m x l*m symmetric matrix
// - b      : l   x n   matrix
// - work   : l*m x n*m matrix
// - result : n*m x n*m symmetric matrix
//
void tensorMatrixCongruenceTranspose(const Matrix &a,
                                     const Matrix &b,
                                     Matrix &work,
                                     Matrix &result) {
  int m = a.rows / b.rows;

  assert(result.rows == b.cols * m);
  assert(result.cols == b.cols * m);

  assert(work.rows == a.rows);
  assert(work.cols == result.cols);

  // work = a b'
  for (int c = 0; c < work.cols; c++) {
    int bCol       = c % b.cols;
    int aColOffset = (c / b.cols) * b.rows;

    for (int r = 0; r < work.rows; r++) {

      Real tmp = 0;
      for (int k = 0; k < b.rows; k++) {
        tmp += a.elt(r, aColOffset + k) * b.elt(k, bCol);
      }

      work.elt(r, c) = tmp;
    }
  }

  // result = b'^T work
  for (int c = 0; c < result.cols; c++) {

    // since result is symmetric, only compute its upper triangle
    for (int r = 0; r <= c; r++) {
      int bCol          = r % b.cols;
      int workRowOffset = (r / b.cols) * b.rows;

      Real tmp = 0;
      for (int k = 0; k < b.rows; k++) {
        tmp += b.elt(k, bCol) * work.elt(workRowOffset + k, c);
      }

      result.elt(r, c) = tmp;

      // lower triangle is the same as upper triangle
      if (c != r) {
        result.elt(c, r) = tmp;
      }
    }
  }
}

// result = B'^T A^{-1} B', where B' = B \otimes 1
// Inputs:
// - L      : l*m x l*m cholesky decomposition of A
// - B      : l   x n   matrix
// - Work   : l*m x n*m matrix
// - Result : n*m x n*m symmetric matrix
//
void tensorMatrixInvCongruenceTransposeWithCholesky(const Matrix &L,
                                                    const Matrix &B,
                                                    Matrix &Work,
                                                    Matrix &Result) {
  // Work = L^{-1} (B \otimes 1);
  for (int cw = 0; cw < Work.cols; cw++) {
    int mc  = cw / B.cols;

    for (int rw = mc*B.rows; rw < Work.rows; rw++) {
      int mr = rw / B.rows;

      Real tmp = (mr != mc) ? Real(0) : B.elt(rw % B.rows, cw % B.cols);
      for (int cl = mc*B.rows; cl < rw; cl++)
        tmp -= L.elt(rw, cl)*Work.elt(cl, cw);

      Work.elt(rw, cw) = tmp/L.elt(rw, rw);
    }
  }

  // Result = Work^T Work
  for (int cr = 0; cr < Result.cols; cr++) {
    int mc = cr / B.cols;

    for (int rr = 0; rr <= cr; rr++) {
      int mr = rr / B.cols;

      Real tmp = 0;
      for (int rw = max(mr, mc)*B.rows; rw < Work.rows; rw++)
        tmp += Work.elt(rw, cr)*Work.elt(rw, rr);

      Result.elt(rr, cr) = tmp;
      if (rr != cr)
        Result.elt(cr, rr) = tmp;
    }
  }
}

void computeBilinearPairings(const BlockDiagonalMatrix &A,
                             const vector<Matrix> &bilinearBases,
                             vector<Matrix> &workspace,
                             BlockDiagonalMatrix &result) {
  #pragma omp parallel for schedule(dynamic)
  for (unsigned int b = 0; b < bilinearBases.size(); b++)
    tensorMatrixCongruenceTranspose(A.blocks[b], bilinearBases[b], workspace[b], result.blocks[b]);
}

void computeInvBilinearPairingsWithCholesky(const BlockDiagonalMatrix &L,
                                            const vector<Matrix> &bilinearBases,
                                            vector<Matrix> &workspace,
                                            BlockDiagonalMatrix &result) {
  #pragma omp parallel for schedule(dynamic)
  for (unsigned int b = 0; b < bilinearBases.size(); b++)
    tensorMatrixInvCongruenceTransposeWithCholesky(L.blocks[b], bilinearBases[b], workspace[b], result.blocks[b]);
}

// result = V D V^T, where D=diag(d) is a diagonal matrix
// Inputs:
// - d        : pointer to beginning of a length-V.cols vector
// - V        : V.rows x V.cols Matrix
// - blockRow : integer < k
// - blockCol : integer < k
// - result   : (k*V.rows) x (k*V.rows) square Matrix
//
void diagonalCongruence(Real const *d,
                        const Matrix &V,
                        const int blockRow,
                        const int blockCol,
                        Matrix &result) {

  for (int p = 0; p < V.rows; p++) {
    for (int q = 0; q <= p; q++) {
      Real tmp = 0;

      for (int n = 0; n < V.cols; n++)
        tmp += *(d+n) * V.elt(p, n)*V.elt(q, n);
      
      result.elt(blockRow*V.rows + p, blockCol*V.rows + q) = tmp;
      if (p != q)
        result.elt(blockRow*V.rows + q, blockCol*V.rows + p) = tmp;
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
        for (vector<int>::const_iterator b = sdp.blocks[j].begin(); b != sdp.blocks[j].end(); b++) {
          tmp += (BilinearPairingsXInv.blocks[*b].elt(ej_s1 + k1, ej_r2 + k2)*BilinearPairingsY.blocks[*b].elt(ej_s2 + k2, ej_r1 + k1) +
                  BilinearPairingsXInv.blocks[*b].elt(ej_r1 + k1, ej_r2 + k2)*BilinearPairingsY.blocks[*b].elt(ej_s2 + k2, ej_s1 + k1) +
                  BilinearPairingsXInv.blocks[*b].elt(ej_s1 + k1, ej_s2 + k2)*BilinearPairingsY.blocks[*b].elt(ej_r2 + k2, ej_r1 + k1) +
                  BilinearPairingsXInv.blocks[*b].elt(ej_r1 + k1, ej_s2 + k2)*BilinearPairingsY.blocks[*b].elt(ej_r2 + k2, ej_s1 + k1))/4;
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
      for (vector<int>::const_iterator b = sdp.blocks[j].begin(); b != sdp.blocks[j].end(); b++) {
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

void constraintMatrixWeightedSum(const SDP &sdp, const Vector x, BlockDiagonalMatrix &result)  {
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

      for (vector<int>::const_iterator b = sdp.blocks[j].begin(); b != sdp.blocks[j].end(); b++) {
        diagonalCongruence(&x[p], sdp.bilinearBases[*b], r, s, result.blocks[*b]);
        
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
      for (vector<int>::const_iterator b = sdp.blocks[j].begin(); b != sdp.blocks[j].end(); b++) {

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
  timers["schurblocks/cholesky"].resume();
  computeSchurBlocks(sdp, BilinearPairingsXInv, BilinearPairingsY, SchurBlocks);
  choleskyDecompositionStabilized(SchurBlocks, SchurBlocksCholesky,
                                  schurStabilizeIndices,
                                  schurStabilizeLambdas,
                                  cast2double(choleskyStabilizeThreshold));
  timers["schurblocks/cholesky"].stop();
  // SchurUpdateLowRank = {{- 1, 0}, {E, G}}
  timers["make schur update"].resume();
  SchurUpdateLowRank.copyFrom(sdp.FreeVarMatrix);  
  blockMatrixLowerTriangularSolve(SchurBlocksCholesky, SchurUpdateLowRank);
  int updateColumns = SchurUpdateLowRank.cols;

  stabilizeBlockIndices.resize(0);
  stabilizeBlockUpdateRow.resize(0);
  stabilizeBlockUpdateColumn.resize(0);

  for (unsigned int b = 0; b < SchurBlocks.blocks.size(); b++) {
    if (schurStabilizeIndices[b].size() > 0) {
      int startIndex = schurStabilizeIndices[b][0];
      int blockRows  = SchurBlocks.blocks[b].rows - startIndex;
      int blockCols  = schurStabilizeIndices[b].size();

      stabilizeBlockIndices.push_back(b);
      stabilizeBlockUpdateRow.push_back(SchurBlocks.blockStartIndices[b] + startIndex);
      stabilizeBlockUpdateColumn.push_back(updateColumns);
      updateColumns += blockCols;

      stabilizeBlocks[b].setRowsCols(blockRows, blockCols);
      stabilizeBlocks[b].setZero();
      for (unsigned int c = 0; c < schurStabilizeIndices[b].size(); c++) {
        int r = schurStabilizeIndices[b][c] - startIndex;
        stabilizeBlocks[b].elt(r, c) = schurStabilizeLambdas[b][c];
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
  timers["make schur update"].stop();
  timers["make Q"].resume();
  // Q = SchurUpdateLowRank^T SchurUpdateLowRank - {{0,0},{0,1}}
  Q.setRowsCols(updateColumns, updateColumns);
  Q.setZero();

  matrixSquareIntoBlock(SchurUpdateLowRank, Q, 0, 0);

  // LowerRight(Q) = G^T G - 1
  for (unsigned int j = 0; j < stabilizeBlockIndices.size(); j++) {
    int b = stabilizeBlockIndices[j];
    int c = stabilizeBlockUpdateColumn[j];
    matrixSquareIntoBlock(stabilizeBlocks[b], Q, c, c);
    for (int i = c; i < c + stabilizeBlocks[b].cols; i++)
      Q.elt(i,i) -= 1;
  }

  // LowerLeft(Q) = G^T U
  # pragma omp parallel for schedule(dynamic)
  for (unsigned int j = 0; j < stabilizeBlockIndices.size(); j++) {
    int b = stabilizeBlockIndices[j];
    int p = stabilizeBlockUpdateRow[j];
    int r = stabilizeBlockUpdateColumn[j];
    Rgemm("Transpose", "NoTranspose",
          stabilizeBlocks[b].cols,
          SchurUpdateLowRank.cols,
          stabilizeBlocks[b].rows,
          1,
          &stabilizeBlocks[b].elements[0],
          stabilizeBlocks[b].rows,
          &SchurUpdateLowRank.elt(p, 0),
          SchurUpdateLowRank.rows,
          0,
          &Q.elt(r, 0),
          Q.rows);
  }

  // UpperRight(Q) = LowerLeft(Q)^T
  # pragma omp parallel for schedule(static)
  for (int c = 0; c < SchurUpdateLowRank.cols; c++)
    for (int r = SchurUpdateLowRank.cols; r < Q.rows; r++)
      Q.elt(c,r) = Q.elt(r,c);
  timers["make Q"].stop();

  timers["LU of Q"].resume();
  Qpivots.resize(Q.rows);
  LUDecomposition(Q, Qpivots);
  timers["LU of Q"].stop();
}


// As inputs, dx and dy are the residues r_x and r_y on the right-hand
// side of the Schur complement equation. As outputs, they are the
// values for dx and dy.
//
void SDPSolver::solveSchurComplementEquation(Vector &dx, Vector &dy) {
  // dx = SchurBlocksCholesky^{-1} dx
  blockMatrixLowerTriangularSolve(SchurBlocksCholesky, dx);

  basicKernelCoords.resize(Q.rows);
  // k_1 = r_y - SchurUpdateLowRank^T dx
  for (unsigned int n = 0; n < dy.size(); n++)
    basicKernelCoords[n] = dy[n];

  vectorScaleMatrixMultiplyTransposeAdd(-1, SchurUpdateLowRank, dx, 1, basicKernelCoords);
  
  // k_2 = -G^T dx
  for (unsigned int j = 0; j < stabilizeBlockIndices.size(); j++) {
    int b = stabilizeBlockIndices[j];
    int pTopLeft = stabilizeBlockUpdateRow[j];
    int cTopLeft = stabilizeBlockUpdateColumn[j];

    for (int c = 0; c < stabilizeBlocks[b].cols; c++) {
      basicKernelCoords[cTopLeft + c] = 0;
      for (int r = 0; r < stabilizeBlocks[b].rows; r++)
        basicKernelCoords[cTopLeft + c] -= dx[pTopLeft + r] * stabilizeBlocks[b].elt(r,c);
    }
  }

  // k = Q^{-1} k
  solveWithLUDecomposition(Q, Qpivots, basicKernelCoords);

  // dx = dx + SchurUpdateLowRank k_1
  vectorScaleMatrixMultiplyAdd(1, SchurUpdateLowRank, basicKernelCoords, 1, dx);
  // dx += G k_2
  for (unsigned int j = 0; j < stabilizeBlockIndices.size(); j++) {
    int b = stabilizeBlockIndices[j];
    int pTopLeft = stabilizeBlockUpdateRow[j];
    int cTopLeft = stabilizeBlockUpdateColumn[j];

    for (int c = 0; c < stabilizeBlocks[b].cols; c++)
      for (int r = 0; r < stabilizeBlocks[b].rows; r++)
        dx[pTopLeft + r] += basicKernelCoords[cTopLeft + c] * stabilizeBlocks[b].elt(r,c);
  }

  // dx = SchurBlocksCholesky^{-T} dx
  blockMatrixLowerTriangularTransposeSolve(SchurBlocksCholesky, dx);
  // dy = k_1
  for (unsigned int n = 0; n < dy.size(); n++)
    dy[n] = basicKernelCoords[n];
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

void SDPSolver::initialize(const SDPSolverParameters &parameters) {
  fillVector(x, 0);
  X.setZero();
  X.addDiagonal(parameters.initialMatrixScalePrimal);
  fillVector(y, 0);
  Y.setZero();
  Y.addDiagonal(parameters.initialMatrixScaleDual);
}

SDPSolverTerminateReason SDPSolver::run(const SDPSolverParameters &parameters,
                                        const path checkpointFile) {
  Real primalStepLength;
  Real dualStepLength;
  
  for (int iteration = 1;; iteration++) {

    if (timers["Last checkpoint"].elapsed().wall >= parameters.checkpointInterval * 1000000000LL)
      saveCheckpoint(checkpointFile);
    if (timers["Solver runtime"].elapsed().wall >= parameters.maxRuntime * 1000000000LL)
      return MaxRuntimeExceeded;

    timers["cholesky"].resume();
    choleskyDecomposition(X, XCholesky);
    choleskyDecomposition(Y, YCholesky);
    timers["cholesky"].stop();
    timers["bilinear pairings"].resume();
    computeInvBilinearPairingsWithCholesky(XCholesky, sdp.bilinearBases, bilinearPairingsWorkspace, BilinearPairingsXInv);
    computeBilinearPairings(Y, sdp.bilinearBases, bilinearPairingsWorkspace, BilinearPairingsY);
    timers["bilinear pairings"].stop();
    // d_k = c_k - Tr(F_k Y) - (D y)_k
    computeDualResidues(sdp, y, BilinearPairingsY, dualResidues);

    // PrimalResidues = sum_p F_p x_p - X - F_0 (F_0 is zero for now)
    computePrimalResidues(sdp, x, X, PrimalResidues);

    status.primalError     = PrimalResidues.maxAbs();
    status.dualError       = maxAbsVector(dualResidues);
    status.primalObjective = sdp.objectiveConst + dotProduct(sdp.primalObjective, x);
    status.dualObjective   = sdp.objectiveConst + dotProduct(sdp.dualObjective, y);

    const bool isPrimalFeasible = status.primalError  < parameters.primalErrorThreshold;
    const bool isDualFeasible   = status.dualError    < parameters.dualErrorThreshold;
    const bool isOptimal        = status.dualityGap() < parameters.dualityGapThreshold;

    if (isPrimalFeasible && isDualFeasible && isOptimal)
      return PrimalDualOptimal;
    else if (isDualFeasible && parameters.findDualFeasible)
      return DualFeasible;
    else if (dualStepLength == 1 && parameters.detectDualFeasibleJump)
      return DualFeasibleJumpDetected;
    // Detect max iterations after computing errors and objective
    // functions for the current point
    else if (iteration > parameters.maxIterations)
      return MaxIterationsExceeded;
    timers["initialize schur solver"].resume();
    initializeSchurComplementSolver(BilinearPairingsXInv, BilinearPairingsY, parameters.choleskyStabilizeThreshold);
    timers["initialize schur solver"].stop();

    Real mu = frobeniusProductSymmetric(X, Y)/X.dim;
    if (mu > parameters.maxComplementarity)
      return MaxComplementarityExceeded;
    timers["predictor"].resume();
    // Mehrotra predictor solution for (dx, dX, dY)
    Real betaPredictor = predictorCenteringParameter(parameters, isPrimalFeasible && isDualFeasible);
    computeSearchDirection(betaPredictor, mu, false);
    timers["predictor"].stop();
    timers["corrector"].resume();
    // Mehrotra corrector solution for (dx, dX, dY)
    Real betaCorrector = correctorCenteringParameter(parameters, X, dX, Y, dY, mu, isPrimalFeasible && isDualFeasible);
    computeSearchDirection(betaCorrector, mu, true);
    timers["corrector"].stop();
    timers["step length"].resume();
    // Step length to preserve positive definiteness
    primalStepLength = stepLength(XCholesky, dX, StepMatrixWorkspace,
                                  eigenvaluesWorkspace, QRWorkspace, parameters);
    dualStepLength   = stepLength(YCholesky, dY, StepMatrixWorkspace,
                                  eigenvaluesWorkspace, QRWorkspace, parameters);
    timers["step length"].stop();
    if (isPrimalFeasible && isDualFeasible) {
      primalStepLength = min(primalStepLength, dualStepLength);
      dualStepLength = primalStepLength;
    }

    printSolverInfo(iteration, mu, status, primalStepLength, dualStepLength, betaCorrector, sdp.dualObjective.size(), Q.rows);

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
