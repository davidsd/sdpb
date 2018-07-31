//=======================================================================
// Copyright 2014-2015 David Simmons-Duffin.
// Distributed under the MIT License.
// (See accompanying file LICENSE or copy at
//  http://opensource.org/licenses/MIT)
//=======================================================================

#pragma once

#include "Block_Diagonal_Matrix.hxx"
#include "Block_Matrix.hxx"
#include "Block_Vector.hxx"
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
  Block_Vector x;

  // a Block_Diagonal_Matrix with block sizes given by
  // sdp.psdMatrixBlockDims()
  Block_Diagonal_Matrix X;

  // a Vector of length N = sdp.dualObjective.size()
  Block_Vector y;

  // a Block_Diagonal_Matrix with the same structure as X
  Block_Diagonal_Matrix Y;

  /********************************************/
  // Solver status

  // NB: here, primalObjective and dualObjective refer to the current
  // values of the objective functions.  In the class SDP, they refer
  // to the vectors c and b.  Hopefully the name-clash won't cause
  // confusion.
  El::BigFloat primal_objective, // f + c . x
    dual_objective,              // f + b . y
    duality_gap;                 // normalized difference of objectives

  // Discrepancy in the primal equality constraints, a
  // Block_Diagonal_Matrix with the same structure as X, called 'P' in
  // the manual:
  //
  //   PrimalResidues = \sum_p A_p x_p - X
  //
  Block_Diagonal_Matrix primal_residues;
  El::BigFloat primal_error; // maxAbs(PrimalResidues)

  // Discrepancy in the dual equality constraints, a Vector of length
  // P, called 'd' in the manual:
  //
  //   dualResidues = c - Tr(A_* Y) - B y
  //
  Block_Vector dual_residues;
  El::BigFloat dual_error; // maxAbs(dualResidues)

  // Create a new solver for a given SDP, with the given parameters
  SDP_Solver(const boost::filesystem::path &sdp_directory,
             const SDP_Solver_Parameters &parameters);

  // Run the solver, backing up to checkpointFile
  SDP_Solver_Terminate_Reason
  run(const boost::filesystem::path checkpoint_file, const bool &debug);

  // Input/output
  void save_checkpoint(const boost::filesystem::path &checkpoint_file);
  void load_checkpoint(const boost::filesystem::path &checkpoint_file);
  void save_solution(const SDP_Solver_Terminate_Reason,
                     const boost::filesystem::path &out_file);
  void print_header();

  void print_iteration(const int &iteration, El::BigFloat &mu,
                       const El::BigFloat &primal_step_length,
                       const El::BigFloat &dual_step_length,
                       const El::BigFloat &beta_corrector);

  void
  test_multiplication(const int m_init, const int m_fin, const int m_step);

private:
  // Compute data needed to solve the Schur complement equation
  void initialize_schur_complement_solver(
    const Block_Diagonal_Matrix &bilinear_pairings_X_inv,
    const Block_Diagonal_Matrix &bilinear_pairings_Y,
    const std::vector<size_t> &block_sizes,
    const std::vector<size_t> &block_indices, const size_t &num_schur_blocks,
    const El::Grid &block_grid,
    const bool &debug, Block_Diagonal_Matrix &schur_complement_cholesky,
    Block_Matrix &schur_off_diagonal, El::DistMatrix<El::BigFloat> &Q);

  // Compute (dx, dX, dy, dY), given the current mu, a reduction
  // parameter beta.  `correctorPhase' specifies whether to use the
  // R-matrix corresponding to the corrector step (if false, we use
  // the predictor R-matrix)
  void compute_search_direction(
    const Block_Diagonal_Matrix &schur_complement_cholesky,
    const Block_Matrix &schur_off_diagonal,
    const Block_Diagonal_Matrix &X_cholesky, const El::BigFloat beta,
    const El::BigFloat &mu, const bool correctorPhase,
    const El::DistMatrix<El::BigFloat> &Q, Block_Vector &dx,
    Block_Diagonal_Matrix &dX, Block_Vector &dy, Block_Diagonal_Matrix &dY);
};
