//=======================================================================
// Copyright 2014 David Simmons-Duffin.
// Distributed under the MIT License.
// (See accompanying file LICENSE or copy at
//  http://opensource.org/licenses/MIT)
//=======================================================================


#ifndef SDPB_SDPSOLVER_H_
#define SDPB_SDPSOLVER_H_

#include <iostream>
#include <ostream>
#include <vector>
#include "boost/filesystem.hpp"
#include "types.h"
#include "Vector.h"
#include "Matrix.h"
#include "BlockDiagonalMatrix.h"
#include "SDP.h"

using std::vector;
using std::ostream;
using std::endl;
using boost::filesystem::path;

// Parameters to control the behavior of an SDPSolver. See the manual
// for a detailed description of each.
//
class SDPSolverParameters {
public:
  int maxIterations;
  int maxRuntime;
  int checkpointInterval;
  bool noFinalCheckpoint;
  bool findPrimalFeasible;
  bool findDualFeasible;
  bool detectPrimalFeasibleJump;
  bool detectDualFeasibleJump;
  int precision;
  int maxThreads;
  Real dualityGapThreshold;
  Real primalErrorThreshold;
  Real dualErrorThreshold;
  Real initialMatrixScalePrimal;
  Real initialMatrixScaleDual;
  Real feasibleCenteringParameter;
  Real infeasibleCenteringParameter;
  Real stepLengthReduction;
  Real choleskyStabilizeThreshold;
  Real maxComplementarity;

  // Set the precision of all Real parameters to equal 'precision'.
  // This is necessary because 'precision' might be set (via the
  // command line or a file) after initializing other parameters.
  //
  void resetPrecision() {
    setPrecision(dualityGapThreshold,          precision);
    setPrecision(primalErrorThreshold,         precision);
    setPrecision(dualErrorThreshold,           precision);
    setPrecision(initialMatrixScalePrimal,     precision);
    setPrecision(initialMatrixScaleDual,       precision);
    setPrecision(feasibleCenteringParameter,   precision);
    setPrecision(infeasibleCenteringParameter, precision);
    setPrecision(stepLengthReduction,          precision);
    setPrecision(choleskyStabilizeThreshold,   precision);
    setPrecision(maxComplementarity,           precision);
  }

  friend ostream& operator<<(ostream& os, const SDPSolverParameters& p);
};

// Reasons for terminating the solver.  See the manual for a detailed
// description of each.
//
enum SDPSolverTerminateReason {
  PrimalDualOptimal,
  PrimalFeasible,
  DualFeasible,
  PrimalFeasibleJumpDetected,
  DualFeasibleJumpDetected,
  MaxComplementarityExceeded,
  MaxIterationsExceeded,
  MaxRuntimeExceeded,
};

ostream &operator<<(ostream& os, const SDPSolverTerminateReason& r);

// SDPSolverStatus contains the information needed to determine
// whether to terminate or not.
//
class SDPSolverStatus {
public:
  Real primalObjective; // f + c . x
  Real dualObjective;   // f + b . y
  Real primalError;     // maxAbs(PrimalResidues)
  Real dualError;       // maxAbs(dualResidues)

  Real dualityGap() const {
    return abs(primalObjective - dualObjective) /
      max(Real(abs(primalObjective) + abs(dualObjective)), Real(1));
  }

  friend ostream& operator<<(ostream& os, const SDPSolverStatus& s);
};

// SDPSolver contains the data structures needed during the running of
// the interior point algorithm.  Each structure is allocated when an
// SDPSolver is initialized, and reused in each iteration.
//
class SDPSolver {
public:
  // SDP to solve.
  SDP sdp;

  // Objective values and errors, re-evaluated each iteration
  SDPSolverStatus status;

  /********************************************/
  // The current point.

  // a Vector of length P = sdp.primalObjective.size()
  Vector x;

  // a BlockDiagonalMatrix with structure given by sdp.psdMatrixBlockDims()
  BlockDiagonalMatrix X;

  // a Vector of length N = sdp.dualObjective.size()
  Vector y;

  // a BlockDiagonalMatrix with the same structure as X
  BlockDiagonalMatrix Y;

  /********************************************/
  // The search direction.
  //
  // These quantities have the same structure as (x, X, y, Y). They
  // are computed twice each iteration: once in the predictor step,
  // and once in the corrector step.
  //
  Vector dx;
  BlockDiagonalMatrix dX;
  Vector dy;
  BlockDiagonalMatrix dY;

  // Discrepancy in the dual equality constraints, a Vector of length
  // P, called 'd' in the manual:
  //
  //   dualResidues = c - Tr(A_* Y) - B y
  //
  Vector dualResidues;

  // Discrepancy in the primal equality constraints, a
  // BlockDiagonalMatrix with the same structure as X, called 'P' in
  // the manual:
  //
  //   PrimalResidues = \sum_P A_p x_p - X
  //
  BlockDiagonalMatrix PrimalResidues;

  /********************************************/
  // Intermediate computations.

  // the Cholesky decompositions of X and Y, each
  // BlockDiagonalMatrices with the same structure as X and Y
  BlockDiagonalMatrix XCholesky;
  BlockDiagonalMatrix YCholesky;

  // Z = X^{-1} (PrimalResidues Y - R), a BlockDiagonalMatrix with the
  // same structure as X and Y
  BlockDiagonalMatrix Z;

  // R = mu I - X Y for the predictor step
  // R = mu I - X Y - dX dY for the corrector step
  // R is a BlockDiagonalMatrix with the same structure as X and Y.
  BlockDiagonalMatrix R;

  // Bilinear pairings needed for computing the Schur complement
  // matrix.  For example,
  //
  //   BilinearPairingsXInv.blocks[b].elt(
  //     (d_j+1) s + k1,
  //     (d_j+1) r + k2
  //   ) = \chi_{b,k1}^T (X.blocks[b]^{-1})^{(s,r)} \chi_{b,k2}
  //
  //     0 <= k1,k2 <= sdp.degrees[j] = d_j
  //     0 <= s,r < sdp.dimensions[j] = m_j
  //
  // where j corresponds to b and M^{(s,r)} denotes the (s,r)-th
  // (d_j+1)x(d_j+1) block of M.
  //
  // BilinearPairingsXInv has one block for each block of X.  The
  // dimension of BilinearPairingsXInv.block[b] is (d_j+1)*m_j.  See
  // SDP.h for more information on d_j and m_j.
  // 
  BlockDiagonalMatrix BilinearPairingsXInv;
  //
  // BilinearPairingsY is analogous to BilinearPairingsXInv, with
  // X^{-1} -> Y.
  //
  BlockDiagonalMatrix BilinearPairingsY;

  // The Schur complement matrix S: a BlockDiagonalMatrix with one
  // block for each 0 <= j < J.  SchurBlocks.blocks[j] has dimension
  // (d_j+1)*m_j*(m_j+1)/2
  //
  BlockDiagonalMatrix SchurBlocks;

  // SchurBlocksCholesky = L': the Cholesky decomposition of the
  // stabilized Schur complement matrix S' = S + U U^T.
  BlockDiagonalMatrix SchurBlocksCholesky;

  // SchurOffDiagonal = L'^{-1} FreeVarMatrix: needed in solving the
  // Schur complement equation.
  Matrix SchurOffDiagonal;

  // As explained in the manual, we use the `stabilized' Schur
  // complement matrix S' = S + U U^T.  Here, the 'update' matrix U
  // has columns given by
  //
  //   U = ( Lambda_{p_1} e_{p_1}, ..., Lambda_{p_M} e_{p_M} )
  //
  // where e_p is a unit vector in the p-th direction and the
  // Lambda_{p_m} are constants.  If p_m appears above, we say the
  // direction p_m has been `stabilized.'  Because U is sparse, we
  // encode it in the smaller data structures below.
  //
  // Recall: S is block diagonal, with blocks labeled by 0 <= j < J.
  //
  // schurStabilizeIndices[j] = a list of which directions of
  // SchurBlocks.blocks[j] have been stabilized (counting from 0 at
  // the upper left of each block), for 0 <= j < J.
  vector<vector<int> > schurStabilizeIndices;
  //
  // schurStabilizeLambdas[j] = a list of constants Lambda for each
  // stabilized direction in schurStabilizeIndices[j], for 0 <= j < J.
  vector<vector<Real> > schurStabilizeLambdas;
  //
  // a list of block indices {j_0, j_1, ..., j_{M-1} } for blocks
  // which have at least one stabilized direction.  We say the blocks
  // j_m have been `stabilized.'  Note that the number M of stabilized
  // blocks could change with each iteration.
  vector<int> stabilizeBlockIndices;
  //
  // For each 0 <= m < M, stabilizeBlockUpdateRow records the row of
  // U corresponding to the first stabilized direction in the j_m'th
  // block of SchurBlocks:
  //
  //   stabilizeBlockUpdateRow[m] = schurStabilizeIndices[j_m][0] +
  //                                SchurBlocks.blockStartIndices[j_m]
  // 
  vector<int> stabilizeBlockUpdateRow;
  //
  // For each 0 <= m < M, stabilizeBlockUpdateColumn records the
  // column of U corresponding to the first stabilized direction in
  // the j_m'th block of SchurBlocks.
  vector<int> stabilizeBlockUpdateColumn;
  //
  // For each 0 <= m < M, stabilizeBlocks is the submatrix of U
  // corresponding to the j_m'th block of SchurBlocks.
  vector<Matrix> stabilizeBlocks;

  // Q = B' L'^{-T} L'^{-1} B' - {{0, 0}, {0, 1}}, where B' =
  // (FreeVarMatrix U).  Q is needed in the factorization of the Schur
  // complement equation.  Q has dimension N'xN', where
  //
  //   N' = cols(B) + cols(U) = N + cols(U)
  //
  // where N is the dimension of the dual objective function.  Note
  // that N' could change with each iteration.
  Matrix Q;
  // a vector of length N', needed for the LU decomposition of Q.
  vector<Integer> Qpivots;

  // dyExtended = (dy z), where z are the extra coordinates introduced
  // in the stabilized Schur complement equation.
  Vector dyExtended;

  /********************************************/
  // Additional workspace variables

  // A BlockDiagonalMatrix with the same structure as X, Y.  Needed
  // for computing step lengths which preserve positive
  // semidefiniteness.
  BlockDiagonalMatrix StepMatrixWorkspace;

  // 
  vector<Matrix> bilinearPairingsWorkspace;
  vector<Vector> eigenvaluesWorkspace;
  vector<Vector> QRWorkspace;

  SDPSolver(const SDP &sdp);
  void initialize(const SDPSolverParameters &parameters);
  SDPSolverTerminateReason run(const SDPSolverParameters &parameters, const path checkpointFile);
  void initializeSchurComplementSolver(const BlockDiagonalMatrix &BilinearPairingsXInv,
                                       const BlockDiagonalMatrix &BilinearPairingsY,
                                       const Real &choleskyStabilizeThreshold);
  void solveSchurComplementEquation(Vector &dx, Vector &dz);
  void computeSearchDirection(const Real &beta, const Real &mu, const bool correctorPhase);
  void saveCheckpoint(const path &checkpointFile);
  void loadCheckpoint(const path &checkpointFile);
  void saveSolution(const SDPSolverTerminateReason, const path &outFile);
};

void printSolverHeader();
void printSolverInfo(int iteration,
                     Real mu,
                     SDPSolverStatus status,
                     Real primalStepLength,
                     Real dualStepLength,
                     Real betaCorrector,
                     int dualObjectiveSize,
                     int Qrows);

#endif  // SDPB_SDPSOLVER_H_
