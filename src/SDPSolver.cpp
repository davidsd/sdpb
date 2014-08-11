#include <iostream>
#include <ostream>
#include "omp.h"
#include "boost/archive/text_oarchive.hpp"
#include "boost/archive/text_iarchive.hpp"
#include "boost/filesystem.hpp"
#include "boost/filesystem/fstream.hpp"
#include "SDPSolver.h"
#include "serialize.h"
#include "Timers.h"

using boost::filesystem::path;
using boost::timer::nanosecond_type;
using std::cout;

SDPSolver::SDPSolver(const SDP &sdp):
  sdp(sdp),
  x(sdp.primalObjective.size(), 0),
  X(sdp.psdMatrixBlockDims()),
  Y(X),
  dx(x),
  dX(X),
  dY(Y),
  dualResidues(x),
  dualResiduesReduced(sdp.primalObjective.size() - sdp.dualObjective.size()),
  PrimalResidues(X),
  FreeVarMatrixReduced(sdp.primalObjective.size() - sdp.dualObjective.size(), sdp.dualObjective.size()),
  dualObjectiveReduced(sdp.dualObjective.size()),
  XCholesky(X),
  YCholesky(X),
  Z(X),
  R(X),
  BilinearPairingsXInv(sdp.bilinearPairingBlockDims()),
  BilinearPairingsY(BilinearPairingsXInv),
  SchurBlocks(sdp.schurBlockDims()),
  SchurBlocksCholesky(SchurBlocks),
  SchurUpdateLowRank(sdp.FreeVarMatrix),
  Q(sdp.FreeVarMatrix.cols, sdp.FreeVarMatrix.cols),
  Qpivots(sdp.FreeVarMatrix.cols),
  basicKernelCoords(Q.rows),
  BasicKernelSpan(sdp.FreeVarMatrix),
  schurStabilizeIndices(SchurBlocks.blocks.size()),
  schurStabilizeLambdas(SchurBlocks.blocks.size()),
  schurStabilizeVectors(SchurBlocks.blocks.size()),
  StepMatrixWorkspace(X)
{
  // initialize bilinearPairingsWorkspace, eigenvaluesWorkspace, QRWorkspace 
  for (unsigned int b = 0; b < sdp.bilinearBases.size(); b++) {
    bilinearPairingsWorkspace.push_back(Matrix(X.blocks[b].rows, BilinearPairingsXInv.blocks[b].cols));
    eigenvaluesWorkspace.push_back(Vector(X.blocks[b].rows));
    QRWorkspace.push_back(Vector(3*X.blocks[b].rows - 1));
  }

  basicIndices = linearlyIndependentRowIndices(sdp.FreeVarMatrix);
  for (int i = 0, p = 0; p < sdp.FreeVarMatrix.rows; p++)
    if (p == basicIndices[i])
      i++;
    else
      nonBasicIndices.push_back(p);

  // Computations needed for free variable elimination
  Matrix DBLU(sdp.dualObjective.size(), sdp.dualObjective.size());
  vector<Integer> DBLUpivots(sdp.dualObjective.size());

  // LU Decomposition of D_B
  for (int n = 0; n < DBLU.cols; n++)
    for (int m = 0; m < DBLU.rows; m++)
      DBLU.elt(m,n) = sdp.FreeVarMatrix.elt(basicIndices[m],n);
  LUDecomposition(DBLU, DBLUpivots);

  // Compute E = - D_N D_B^{-1}
  // ET = -D_N^T
  Matrix FreeVarMatrixReducedT(FreeVarMatrixReduced.cols, FreeVarMatrixReduced.rows);
  for (int p = 0; p < FreeVarMatrixReducedT.cols; p++)
    for (int n = 0; n < FreeVarMatrixReducedT.rows; n++)
      FreeVarMatrixReducedT.elt(n, p) = -sdp.FreeVarMatrix.elt(nonBasicIndices[p], n);
  // ET = D_B^{-1 T} ET = -D_B^{-1 T} D_N^T
  solveWithLUDecompositionTranspose(DBLU, DBLUpivots,
                                    &FreeVarMatrixReducedT.elements[0],
                                    FreeVarMatrixReducedT.cols,
                                    FreeVarMatrixReducedT.rows);
  // E = ET^T
  transpose(FreeVarMatrixReducedT, FreeVarMatrixReduced);

  // dualObjectiveReduced = D_B^{-T} f
  for (unsigned int n = 0; n < dualObjectiveReduced.size(); n++)
    dualObjectiveReduced[n] = sdp.dualObjective[n];
  solveWithLUDecompositionTranspose(DBLU, DBLUpivots, &dualObjectiveReduced[0], 1, dualObjectiveReduced.size());

  // BasicKernelSpan = ( -1 \\ E)
  BasicKernelSpan.setZero();
  for (int c = 0; c < FreeVarMatrixReduced.cols; c++)
    for (int r = 0; r < FreeVarMatrixReduced.rows; r++)
      BasicKernelSpan.elt(nonBasicIndices[r], c) = FreeVarMatrixReduced.elt(r, c);
  for (int c = 0; c < FreeVarMatrixReduced.cols; c++)
    BasicKernelSpan.elt(basicIndices[c], c) = -1;

  for (unsigned int b = 0; b < SchurBlocks.blocks.size(); b++)
    schurStabilizeVectors[b].resize(SchurBlocks.blocks[b].rows);
}

void printSolverHeader() {
  cout << "\n     mu       P-obj       D-obj     gap         P-err        D-err       P-step   D-step   beta\n";
  cout << "---------------------------------------------------------------------------------------------------\n";
}

void printSolverInfo(int iteration,
                     Real mu,
                     SDPSolverStatus status,
                     bool isPrimalFeasible,
                     bool isDualFeasible,
                     Real primalStepLength,
                     Real dualStepLength,
                     Real betaCorrector) {
  gmp_fprintf(stdout,
              "%3d  %4.1Fe  %+7.2Fe  %+7.2Fe  %+7.2Fe  %s%+7.2Fe%s  %s%+7.2Fe%s  %4.1Fe  %4.1Fe  %4.2Fe\n",
              iteration,
              mu.get_mpf_t(),
              status.primalObjective.get_mpf_t(),
              status.dualObjective.get_mpf_t(),
              status.dualityGap().get_mpf_t(),
              isPrimalFeasible ? "|" : " ", status.primalError.get_mpf_t(), isPrimalFeasible ? "|" : " ",
              isDualFeasible   ? "|" : " ", status.dualError.get_mpf_t(),   isDualFeasible   ? "|" : " ",
              primalStepLength.get_mpf_t(),
              dualStepLength.get_mpf_t(),
              betaCorrector.get_mpf_t());
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
  // X = L^{-1} (B \otimes 1);
  for (int cw = 0; cw < Work.cols; cw++) {
    int mc  = cw / B.cols;

    for (int rw = mc*B.rows; rw < Work.rows; rw++) {
      int mr = rw / B.cols;

      Real tmp = (mr != mc) ? 0 : B.elt(rw % B.rows, cw % B.cols);
      for (int cl = 0; cl < rw; cl++)
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

      for (unsigned int u2 = 0; u2 < sdp.constraintIndices[j].size(); u2++) {
        const int ej_r2 = sdp.constraintIndices[j][u2].r * ej;
        const int ej_s2 = sdp.constraintIndices[j][u2].s * ej;
        const int k2    = sdp.constraintIndices[j][u2].k;

        Real tmp = 0;
        for (vector<int>::const_iterator b = sdp.blocks[j].begin(); b != sdp.blocks[j].end(); b++) {
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

// x_B = g + E^T x_N
void basicCompletion(const Vector &dualObjectiveReduced,
                     const Matrix &FreeVarMatrixReduced,
                     const vector<int> &basicIndices,
                     const vector<int> &nonBasicIndices,
                     Vector &x) {
  assert((int)basicIndices.size()    == FreeVarMatrixReduced.cols);
  assert((int)nonBasicIndices.size() == FreeVarMatrixReduced.rows);
  assert((int)x.size()               == FreeVarMatrixReduced.cols + FreeVarMatrixReduced.rows);

  #pragma omp parallel for schedule(static)
  for (unsigned int n = 0; n < basicIndices.size(); n++) {
    x[basicIndices[n]] = dualObjectiveReduced[n];
    for (unsigned int p = 0; p < nonBasicIndices.size(); p++)
      x[basicIndices[n]] += FreeVarMatrixReduced.elt(p, n) * x[nonBasicIndices[p]];
  }
}

// xReduced_N = x_N + E x_B
void nonBasicShift(const Matrix &FreeVarMatrixReduced,
                   const vector<int> &basicIndices,
                   const vector<int> &nonBasicIndices,
                   const Vector &x,
                   Vector &xReduced) {
  assert((int)basicIndices.size()    == FreeVarMatrixReduced.cols);
  assert((int)nonBasicIndices.size() == FreeVarMatrixReduced.rows);
  assert((int)x.size()               == FreeVarMatrixReduced.cols + FreeVarMatrixReduced.rows);
  assert(nonBasicIndices.size()      == xReduced.size());
  
  #pragma omp parallel for schedule(static)
  for (unsigned int p = 0; p < nonBasicIndices.size(); p++) {
    xReduced[p] = x[nonBasicIndices[p]];
    for (unsigned int n = 0; n < basicIndices.size(); n++)
      xReduced[p] += FreeVarMatrixReduced.elt(p,n) * x[basicIndices[n]];
  }
}

void computeDualResidues(const SDP &sdp,
                         const BlockDiagonalMatrix &Y,
                         const BlockDiagonalMatrix &BilinearPairingsY,
                         Vector &dualResidues) {
  #pragma omp parallel for schedule(dynamic)
  for (unsigned int j = 0; j < sdp.dimensions.size(); j++) {
    const int ej = sdp.degrees[j] +1;

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
      dualResidues[p] += sdp.primalObjective[p];
    }
  }
}

void constraintMatrixWeightedSum(const SDP &sdp, const Vector x, BlockDiagonalMatrix &result)  {
  #pragma omp parallel for schedule(dynamic)
  for (unsigned int j = 0; j < sdp.dimensions.size(); j++) {
    const int dj = sdp.degrees[j];
    int p  = sdp.constraintIndices[j][0].p;

    for (int s = 0; s < sdp.dimensions[j]; s++) {
      for (int r = 0; r <= s; r++) {
        for (vector<int>::const_iterator b = sdp.blocks[j].begin(); b != sdp.blocks[j].end(); b++)
          diagonalCongruence(&x[p], sdp.bilinearBases[*b], r, s, result.blocks[*b]);
        p += dj + 1;
      }
    }
  }

  result.symmetrize();
}

void computeSchurRHS(const SDP &sdp,
                     Vector &dualResidues,
                     BlockDiagonalMatrix &Z, 
                     Vector &r) {

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

Real primalObjectiveValue(const SDP &sdp, const Vector &x) {
  return sdp.objectiveConst + dotProduct(sdp.primalObjective, x);
}

Real dualObjectiveValue(const SDP &sdp, const Vector &dualObjectiveReduced,
                        const vector<int> &basicIndices, const Vector &dualResidues) {
  Real tmp = sdp.objectiveConst;
  for (unsigned int i = 0; i < dualObjectiveReduced.size(); i++)
    tmp += dualObjectiveReduced[i]*dualResidues[basicIndices[i]];
  return tmp;
}

// Implements SDPA's DirectionParameter::MehrotraPredictor
Real predictorCenteringParameter(const SDPSolverParameters &parameters, 
                                 const bool reductionSwitch,
                                 const bool isPrimalDualFeasible) {
  if (isPrimalDualFeasible)
    return 0;
  else if (reductionSwitch)
    return parameters.infeasibleCenteringParameter;
  else
    return 2;
}

// Implements SDPA's DirectionParameter::MehrotraCorrector
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

void addKernelColumn(const Matrix &FreeVarMatrixReduced,
                     const vector<int> &basicIndices,
                     const vector<int> &nonBasicIndices,
                     const int i,
                     const Real &lambda,
                     Matrix &K) {
  K.addColumn();
  int c = K.cols - 1;

  int j = binaryFind(basicIndices.begin(), basicIndices.end(), i) - basicIndices.begin();
  if (j < FreeVarMatrixReduced.cols) {
    for (unsigned int r = 0; r < nonBasicIndices.size(); r++)
      K.elt(nonBasicIndices[r], c) = lambda * FreeVarMatrixReduced.elt(r, j);
  } else {
    K.elt(i, c) = lambda;
  }
}

void SDPSolver::initializeSchurComplementSolver(const BlockDiagonalMatrix &BilinearPairingsXInv,
                                                const BlockDiagonalMatrix &BilinearPairingsY) {

  computeSchurBlocks(sdp, BilinearPairingsXInv, BilinearPairingsY, SchurBlocks);
  choleskyDecomposition(SchurBlocks, SchurBlocksCholesky);

  for (unsigned int b = 0; b < SchurBlocks.blocks.size(); b++)
    stabilizeCholesky(SchurBlocksCholesky.blocks[b],
                      schurStabilizeVectors[b],
                      schurStabilizeIndices[b],
                      schurStabilizeLambdas[b]);
  
  // SchurUpdateLowRank = {{- 1, 0}, {E, G}}
  SchurUpdateLowRank.setCols(BasicKernelSpan.cols);
  SchurUpdateLowRank.copyFrom(BasicKernelSpan);
  for (unsigned int b = 0; b < SchurBlocks.blocks.size(); b++) {
    for (unsigned int i = 0; i < schurStabilizeIndices[b].size(); i++) {
      int fullIndex = SchurBlocks.blockStartIndices[b] + schurStabilizeIndices[b][i];
      addKernelColumn(FreeVarMatrixReduced,
                      basicIndices,
                      nonBasicIndices,
                      fullIndex,
                      schurStabilizeLambdas[b],
                      SchurUpdateLowRank);
    }
    schurStabilizeIndices[b].resize(0);
  }

  // SchurUpdateLowRank = SchurBlocksCholesky^{-1} {{- 1, 0}, {E, G}}
  blockMatrixLowerTriangularSolve(SchurBlocksCholesky, SchurUpdateLowRank);

  // Q = SchurUpdateLowRank^T SchurUpdateLowRank - {{0,0},{0,1}}
  Q.setRowsCols(SchurUpdateLowRank.cols, SchurUpdateLowRank.cols);
  matrixSquare(SchurUpdateLowRank, Q);
  int stabilizerStart = FreeVarMatrixReduced.cols;
  for (int i = stabilizerStart; i < Q.cols; i++)
    Q.elt(i,i) -= 1;

  Qpivots.resize(Q.rows);
  LUDecomposition(Q, Qpivots);
}

void SDPSolver::solveSchurComplementEquation(Vector &dx) {

  // dx = SchurBlocksCholesky^{-1} dx
  blockMatrixLowerTriangularSolve(SchurBlocksCholesky, dx);

  // k = -SchurUpdateLowRank^T dx
  basicKernelCoords.resize(SchurUpdateLowRank.cols);
  vectorScaleMatrixMultiplyTransposeAdd(-1, SchurUpdateLowRank, dx, 0, basicKernelCoords);

  // k = Q^{-1} k
  solveWithLUDecomposition(Q, Qpivots, basicKernelCoords);

  // dx = dx + SchurUpdateLowRank k
  vectorScaleMatrixMultiplyAdd(1, SchurUpdateLowRank, basicKernelCoords, 1, dx);

  // dx = SchurBlocksCholesky^{-T} dx
  blockMatrixLowerTriangularTransposeSolve(SchurBlocksCholesky, dx);
}

void SDPSolver::computeSearchDirection(const Real &beta,
                                       const Real &mu,
                                       const bool correctorPhase) {

  blockDiagonalMatrixScaleMultiplyAdd(-1, X,  Y,  0, R);
  if (correctorPhase)
    blockDiagonalMatrixScaleMultiplyAdd(-1, dX, dY, 1, R);
  R.addDiagonal(beta*mu);

  // Z = Symmetrize(X^{-1} (PrimalResidues Y - R))
  blockDiagonalMatrixMultiply(PrimalResidues, Y, Z);
  Z -= R;
  blockMatrixSolveWithCholesky(XCholesky, Z);
  Z.symmetrize();

  // dx_k = -d_k + Tr(F_k Z)
  computeSchurRHS(sdp, dualResidues, Z, dx);

  // dx_N = B_{NN}^{-1} dx_N, dx_B = E^T dx_N
  solveSchurComplementEquation(dx);

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

SDPSolverTerminateReason SDPSolver::run(const SDPSolverParameters &parameters,
                                        const path checkpointFile) {
  printSolverHeader();
  timers["Run solver"].start();
  timers["Save checkpoint"].start();
  nanosecond_type const checkpointNanoseconds = parameters.checkpointInterval * 1000000000LL;
  nanosecond_type const maxRuntimeNanoseconds = parameters.maxRuntime * 1000000000LL;
  SDPSolverTerminateReason finished = MaxIterationsExceeded;

  for (int iteration = 1;; iteration++) {
    
    if (timers["Save checkpoint"].elapsed().wall >= checkpointNanoseconds) {
      saveCheckpoint(checkpointFile);
      timers["Save checkpoint"].start();
    }
    if (timers["Run solver"].elapsed().wall >= maxRuntimeNanoseconds) {
      finished = MaxRuntimeExceeded;
      break;
    }
    if (iteration > parameters.maxIterations) {
      finished = MaxIterationsExceeded;
      break;
    }

    // Maintain the invariant x_B = g + E^T x_N
    basicCompletion(dualObjectiveReduced, FreeVarMatrixReduced, basicIndices, nonBasicIndices, x);

    choleskyDecomposition(X, XCholesky);
    choleskyDecomposition(Y, YCholesky);

    computeInvBilinearPairingsWithCholesky(XCholesky, sdp.bilinearBases, bilinearPairingsWorkspace, BilinearPairingsXInv);
    computeBilinearPairings(Y, sdp.bilinearBases, bilinearPairingsWorkspace, BilinearPairingsY);

    // d_k = c_k - Tr(F_k Y)
    computeDualResidues(sdp, Y, BilinearPairingsY, dualResidues);
    nonBasicShift(FreeVarMatrixReduced, basicIndices, nonBasicIndices, dualResidues, dualResiduesReduced);

    // PrimalResidues = sum_p F_p x_p - X - F_0 (F_0 is zero for now)
    computePrimalResidues(sdp, x, X, PrimalResidues);

    status.primalError     = PrimalResidues.maxAbs();
    status.dualError       = maxAbsVector(dualResiduesReduced);
    status.primalObjective = primalObjectiveValue(sdp, x);
    status.dualObjective   = dualObjectiveValue(sdp, dualObjectiveReduced, basicIndices, dualResidues);

    const bool isPrimalFeasible = status.primalError  < parameters.primalErrorThreshold;
    const bool isDualFeasible   = status.dualError    < parameters.dualErrorThreshold;
    const bool isOptimal        = status.dualityGap() < parameters.dualityGapThreshold;
    const bool reductionSwitch  = true;

    if (isPrimalFeasible && isDualFeasible && isOptimal) {
      finished = PrimalDualOptimal;
      break;
    } else if (isDualFeasible && status.dualObjective > parameters.maxDualObjective) {
      finished = DualFeasibleMaxObjectiveExceeded;
      break;
    }

    initializeSchurComplementSolver(BilinearPairingsXInv, BilinearPairingsY);

    Real mu = frobeniusProductSymmetric(X, Y)/X.dim;

    // Mehrotra predictor solution for (dx, dX, dY)
    Real betaPredictor = predictorCenteringParameter(parameters, reductionSwitch,
                                                     isPrimalFeasible && isDualFeasible);
    computeSearchDirection(betaPredictor, mu, false);

    // Mehrotra corrector solution for (dx, dX, dY)
    Real betaCorrector = correctorCenteringParameter(parameters, X, dX, Y, dY, mu,
                                                     isPrimalFeasible && isDualFeasible);
    computeSearchDirection(betaCorrector, mu, true);

    // Step length to preserve positive definiteness
    Real primalStepLength = stepLength(XCholesky, dX, StepMatrixWorkspace,
                                       eigenvaluesWorkspace, QRWorkspace, parameters);
    Real dualStepLength   = stepLength(YCholesky, dY, StepMatrixWorkspace,
                                       eigenvaluesWorkspace, QRWorkspace, parameters);

    printSolverInfo(iteration, mu, status, isPrimalFeasible, isDualFeasible,
                    primalStepLength, dualStepLength, betaCorrector);

    // Update current point
    scaleVector(dx, primalStepLength);
    addVector(x, dx);
    dX *= primalStepLength;
    X += dX;
    dY *= dualStepLength;
    Y += dY;
  }
  
  timers["Run solver"].stop();
  saveCheckpoint(checkpointFile);
  timers["Save checkpoint"].start();
  return finished;
}

ostream& operator<<(ostream& os, const SDPSolverParameters& p) {
  os << "maxIterations                = " << p.maxIterations                << endl;
  os << "maxRuntime                   = " << p.maxRuntime                   << endl;
  os << "checkpointInterval           = " << p.checkpointInterval           << endl;
  os << "precision(actual)            = " << p.precision << "(" << mpf_get_default_prec() << ")" << endl;
  os << "maxThreads                   = " << p.maxThreads                   << endl;
  os << "dualityGapThreshold          = " << p.dualityGapThreshold          << endl;
  os << "primalErrorThreshold         = " << p.primalErrorThreshold         << endl;
  os << "dualErrorThreshold           = " << p.dualErrorThreshold           << endl;
  os << "initialMatrixScale           = " << p.initialMatrixScale           << endl;
  os << "feasibleCenteringParameter   = " << p.feasibleCenteringParameter   << endl;
  os << "infeasibleCenteringParameter = " << p.infeasibleCenteringParameter << endl;
  os << "stepLengthReduction          = " << p.stepLengthReduction          << endl;
  os << "maxDualObjective             = " << p.maxDualObjective             << endl;
  return os;
}

ostream &operator<<(ostream& os, const SDPSolverTerminateReason& r) {
  switch(r) {
  case PrimalDualOptimal:
    os << "found primal-dual optimal solution.";
    break;
  case MaxIterationsExceeded:
    os << "maxIterations exceeded.";
    break;
  case MaxRuntimeExceeded:
    os << "maxRuntime exceeded.";
    break;
  case DualFeasibleMaxObjectiveExceeded:
    os << "found dual feasible solution with dualObjective exceeding maxDualObjective.";
    break;
  }
  return os;
}

ostream& operator<<(ostream& os, const SDPSolverStatus& s) {
  os << "primalObjective = " << s.primalObjective << endl;
  os << "dualObjective   = " << s.dualObjective << endl;
  os << "dualityGap      = " << s.dualityGap() << endl;
  os << "primalError     = " << s.primalError << endl;
  os << "dualError       = " << s.dualError << endl;
  return os;
}

void backupCheckpointFile(path const& checkpointFile) {
  path backupFile(checkpointFile);
  backupFile.replace_extension(".ck.bk");
  cout << "Backing up checkpoint to: " << backupFile << endl;
  copy_file(checkpointFile, backupFile, boost::filesystem::copy_option::overwrite_if_exists);
}

void SDPSolver::saveCheckpoint(const path &checkpointFile) {
  if (exists(checkpointFile))
    backupCheckpointFile(checkpointFile);
  boost::filesystem::ofstream ofs(checkpointFile);
  boost::archive::text_oarchive ar(ofs);
  cout << "Saving checkpoint to    : " << checkpointFile << endl;
  boost::serialization::serializeSDPSolverState(ar, x, X, Y);
}

void SDPSolver::loadCheckpoint(const path &checkpointFile) {
  boost::filesystem::ifstream ifs(checkpointFile);
  boost::archive::text_iarchive ar(ifs);
  cout << "Loading checkpoint from : " << checkpointFile << endl;
  boost::serialization::serializeSDPSolverState(ar, x, X, Y);
}

void SDPSolver::initialize(const SDPSolverParameters &parameters) {
  fillVector(x, 0);
  X.setZero();
  X.addDiagonal(parameters.initialMatrixScale);
  Y.setZero();
  Y.addDiagonal(parameters.initialMatrixScale);
}

void SDPSolver::saveSolution(const path &outFile) {
  boost::filesystem::ofstream ofs(outFile);
  cout << "Saving solution to: " << outFile << endl;
  ofs << "foo" << endl;
}
