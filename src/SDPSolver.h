#ifndef SDP_BOOTSTRAP_SDPSOLVER_H_
#define SDP_BOOTSTRAP_SDPSOLVER_H_

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

class SDPSolverParameters {
public:
  int maxIterations;
  int maxRuntime;
  int checkpointInterval;
  int precision;
  int maxThreads;
  Real dualityGapThreshold;
  Real primalErrorThreshold;
  Real dualErrorThreshold;
  Real initialMatrixScale;
  Real feasibleCenteringParameter;
  Real infeasibleCenteringParameter;
  Real stepLengthReduction;
  Real maxDualObjective;

  void resetPrecision() {
    dualityGapThreshold         .set_prec(precision);
    primalErrorThreshold        .set_prec(precision);
    dualErrorThreshold          .set_prec(precision);
    initialMatrixScale          .set_prec(precision);
    feasibleCenteringParameter  .set_prec(precision);
    infeasibleCenteringParameter.set_prec(precision);
    stepLengthReduction         .set_prec(precision);
    maxDualObjective            .set_prec(precision);
  }

  friend ostream& operator<<(ostream& os, const SDPSolverParameters& p);
};

enum SDPSolverTerminateReason {
  PrimalDualOptimal,
  DualFeasibleMaxObjectiveExceeded,
  MaxIterationsExceeded,
  MaxRuntimeExceeded,
};

ostream &operator<<(ostream& os, const SDPSolverTerminateReason& r);

class SDPSolverStatus {
public:
  Real primalObjective;
  Real dualObjective;
  Real primalError;
  Real dualError;

  Real dualityGap() const {
    return abs(primalObjective - dualObjective) /
      max(Real(abs(primalObjective) + abs(dualObjective)), Real(1));
  }

  friend ostream& operator<<(ostream& os, const SDPSolverStatus& s);
};

class SDPSolver {
public:
  SDP sdp;
  SDPSolverStatus status;

  // current point
  Vector x;
  BlockDiagonalMatrix X;
  BlockDiagonalMatrix Y;

  // search direction
  Vector dx;
  BlockDiagonalMatrix dX;
  BlockDiagonalMatrix dY;

  // discrepancies in dual and primal equality constraints
  Vector dualResidues;
  Vector dualResiduesReduced;
  BlockDiagonalMatrix PrimalResidues;

  // For free variable elimination
  Matrix FreeVarMatrixReduced;
  Vector dualObjectiveReduced;
  vector<int> basicIndices;
  vector<int> nonBasicIndices;

  // intermediate computations
  BlockDiagonalMatrix XCholesky;
  BlockDiagonalMatrix YCholesky;
  BlockDiagonalMatrix Z;
  BlockDiagonalMatrix R;
  BlockDiagonalMatrix BilinearPairingsXInv;
  BlockDiagonalMatrix BilinearPairingsY;
  BlockDiagonalMatrix SchurBlocks;
  BlockDiagonalMatrix SchurBlocksCholesky;
  Matrix SchurUpdateLowRank;
  Matrix Q;
  vector<Integer> Qpivots;
  Vector basicKernelCoords;
  Matrix BasicKernelSpan;
  vector<vector<int> > schurStabilizeIndices;
  vector<Real> schurStabilizeLambdas;
  vector<Vector> schurStabilizeVectors;

  // additional workspace variables
  BlockDiagonalMatrix StepMatrixWorkspace;
  vector<Matrix> bilinearPairingsWorkspace;
  vector<Vector> eigenvaluesWorkspace;
  vector<Vector> QRWorkspace;

  SDPSolver(const SDP &sdp);
  void initialize(const SDPSolverParameters &parameters);
  SDPSolverTerminateReason run(const SDPSolverParameters &parameters, const path checkpointFile);
  void initializeSchurComplementSolver(const BlockDiagonalMatrix &BilinearPairingsXInv,
                                       const BlockDiagonalMatrix &BilinearPairingsY);
  void solveSchurComplementEquation(Vector &dx);
  void computeSearchDirection(const Real &beta, const Real &mu, const bool correctorPhase);
  void saveCheckpoint(const path &checkpointFile);
  void loadCheckpoint(const path &checkpointFile);
  void saveSolution(const path &outFile);
};

#endif  // SDP_BOOTSTRAP_SDPSOLVER_H_
