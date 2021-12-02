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

#include "Solver_Parameters.hxx"
#include "../Timers.hxx"

#include <boost/filesystem.hpp>

// SDPSolver contains the data structures needed during the running of
// the interior point algorithm.  Each structure is allocated when an
// SDPSolver is initialized, and reused in each iteration.
//
class SDP_Solver
{
public:
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

  // primal_error is max of both primal_residues and p=(b - B^T x)
  El::BigFloat primal_error_P, primal_error_p; // |P| and |p|
  El::BigFloat primal_error() const
  {
    return std::max(primal_error_P, primal_error_p);
  }

  // Discrepancy in the dual equality constraints, a Vector of length
  // P, called 'd' in the manual:
  //
  //   dualResidues = c - Tr(A_* Y) - B y
  //
  Block_Vector dual_residues;
  El::BigFloat dual_error; // maxAbs(dualResidues)

  int64_t current_generation;
  boost::optional<int64_t> backup_generation;

  SDP_Solver(const Solver_Parameters &parameters,
             const Verbosity &verbosity,
             const bool &require_initial_checkpoint,
             const Block_Info &block_info, const El::Grid &grid,
             const size_t &dual_objective_b_height);

  SDP_Solver_Terminate_Reason
  run(const Solver_Parameters &parameters,
      const Verbosity &verbosity,
      const boost::property_tree::ptree &parameter_properties,
      const Block_Info &block_info, const SDP &sdp, const El::Grid &grid,
      Timers &timers);

  void step(
    const Solver_Parameters &parameters, const std::size_t &total_psd_rows,
    const bool &is_primal_and_dual_feasible, const Block_Info &block_info,
    const SDP &sdp, const El::Grid &grid,
    const Block_Diagonal_Matrix &X_cholesky,
    const Block_Diagonal_Matrix &Y_cholesky,
    const std::array<
      std::vector<std::vector<std::vector<El::DistMatrix<El::BigFloat>>>>, 2>
      &A_X_inv,
    const std::array<
      std::vector<std::vector<std::vector<El::DistMatrix<El::BigFloat>>>>, 2>
      &A_Y,
    const Block_Vector &primal_residue_p, El::BigFloat &mu,
    El::BigFloat &beta_corrector, El::BigFloat &primal_step_length,
    El::BigFloat &dual_step_length, bool &terminate_now, Timers &timers);

  void
  save_checkpoint(const boost::filesystem::path &checkpoint_directory,
                  const Verbosity &verbosity,
                  const boost::property_tree::ptree &parameter_properties);
  bool
  load_checkpoint(const boost::filesystem::path &checkpoint_directory,
                  const Block_Info &block_info, const Verbosity &verbosity,
                  const bool &require_initial_checkpoint);
};
