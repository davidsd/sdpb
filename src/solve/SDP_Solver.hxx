//=======================================================================
// Copyright 2014-2015 David Simmons-Duffin.
// Distributed under the MIT License.
// (See accompanying file LICENSE or copy at
//  http://opensource.org/licenses/MIT)
//=======================================================================

#pragma once

#include "Block_Diagonal_Matrix.hxx"
#include "SDP.hxx"
#include "SDP_Solver_Terminate_Reason.hxx"

#include "../SDP_Solver_Parameters.hxx"

#include <boost/filesystem.hpp>

// SDPSolver contains the data structures needed during the running of
// the interior point algorithm.  Each structure is allocated when an
// SDPSolver is initialized, and reused in each iteration.
//
class SDP_Solver
{
public:
  // SDP to solve.
  SDP sdp;

  // parameters for initialization and iteration
  SDP_Solver_Parameters parameters;

  /********************************************/
  // Current point

  // a Vector of length P = sdp.primalObjective.size()
  Vector x;

  // a Block_Diagonal_Matrix with block sizes given by
  // sdp.psdMatrixBlockDims()
  Block_Diagonal_Matrix X;

  // a Vector of length N = sdp.dualObjective.size()
  Vector y;

  // a Block_Diagonal_Matrix with the same structure as X
  Block_Diagonal_Matrix Y;

  /********************************************/
  // Search direction
  //
  // These quantities have the same structure as (x, X, y, Y). They
  // are computed twice each iteration: once in the predictor step,
  // and once in the corrector step.
  //
  Vector dx;
  Block_Diagonal_Matrix dX;
  Vector dy;
  Block_Diagonal_Matrix dY;

  /********************************************/
  // Solver status

  // NB: here, primalObjective and dualObjective refer to the current
  // values of the objective functions.  In the class SDP, they refer
  // to the vectors c and b.  Hopefully the name-clash won't cause
  // confusion.
  Real primal_objective; // f + c . x
  Real dual_objective;   // f + b . y
  Real duality_gap;      // normalized difference of objectives

  // Discrepancy in the primal equality constraints, a
  // Block_Diagonal_Matrix with the same structure as X, called 'P' in
  // the manual:
  //
  //   PrimalResidues = \sum_p A_p x_p - X
  //
  Block_Diagonal_Matrix primal_residues;
  Real primal_error; // maxAbs(PrimalResidues)

  // Discrepancy in the dual equality constraints, a Vector of length
  // P, called 'd' in the manual:
  //
  //   dualResidues = c - Tr(A_* Y) - B y
  //
  Vector dual_residues;
  Real dual_error; // maxAbs(dualResidues)

  /********************************************/
  // Intermediate computations.

  // Z = X^{-1} (PrimalResidues Y - R), a Block_Diagonal_Matrix with the
  // same block sizes as X and Y
  Block_Diagonal_Matrix Z;

  // R = mu I - X Y for the predictor step
  // R = mu I - X Y - dX dY for the corrector step
  // R has the same block sizes as X and Y.
  Block_Diagonal_Matrix R;

  // The Schur complement matrix S: a Block_Diagonal_Matrix with one
  // block for each 0 <= j < J.  SchurComplement.blocks[j] has dimension
  // (d_j+1)*m_j*(m_j+1)/2
  //
  Block_Diagonal_Matrix schur_complement;

  // SchurComplementCholesky = L', the Cholesky decomposition of the
  // Schur complement matrix S.
  Block_Diagonal_Matrix schur_complement_cholesky;

  // SchurOffDiagonal = L'^{-1} FreeVarMatrix, needed in solving the
  // Schur complement equation.
  Matrix schur_off_diagonal;

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
  std::vector<Integer> Q_pivots;

  /********************************************/
  // Methods

  // Create a new solver for a given SDP, with the given parameters
  SDP_Solver(const std::vector<boost::filesystem::path> &sdp_files,
             const SDP_Solver_Parameters &parameters);

  // Run the solver, backing up to checkpointFile
  SDP_Solver_Terminate_Reason
  run(const boost::filesystem::path checkpoint_file);

  // Input/output
  void save_checkpoint(const boost::filesystem::path &checkpoint_file);
  void load_checkpoint(const boost::filesystem::path &checkpoint_file);
  void save_solution(const SDP_Solver_Terminate_Reason,
                     const boost::filesystem::path &out_file);
  void print_header();
  void print_iteration(int iteration, Real mu, Real primal_step_length,
                       Real dual_step_length, Real beta_corrector);

  void
  test_multiplication(const int m_init, const int m_fin, const int m_step);

private:
  // Compute data needed to solve the Schur complement equation
  void initialize_schur_complement_solver(
    const Block_Diagonal_Matrix &bilinear_pairings_X_Inv,
    const Block_Diagonal_Matrix &bilinear_pairings_Y);

  // Solve the Schur complement equation in-place for dx, dy.  dx and
  // dy will initially be set to the right-hand side of the equation.
  // This function replaces them with the corresponding solutions.
  void solve_schur_complement_equation(Vector &dx, Vector &dy);

  // Compute (dx, dX, dy, dY), given the current mu, a reduction
  // parameter beta.  `correctorPhase' specifies whether to use the
  // R-matrix corresponding to the corrector step (if false, we use
  // the predictor R-matrix)
  void compute_search_direction(const Block_Diagonal_Matrix &X_cholesky,
                                const Real &beta, const Real &mu,
                                const bool correctorPhase);
};
