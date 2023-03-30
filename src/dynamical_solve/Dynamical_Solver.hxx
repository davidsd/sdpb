//=======================================================================
// Copyright 2014-2015 David Simmons-Duffin.
// Distributed under the MIT License.
// (See accompanying file LICENSE or copy at
//  http://opensource.org/licenses/MIT)
//=======================================================================

#pragma once

#include "../sdp_solve/Block_Diagonal_Matrix.hxx"
#include "../sdp_solve/SDP.hxx"
#include "../sdp_solve/Solver_Parameters.hxx"
#include "Dynamical_Solver_Parameters.hxx"
#include "Dynamical_Solver_Terminate_Reason.hxx"

#include "../Timers.hxx"

#include <boost/filesystem.hpp>

// Dynamical Solver contains the data structures needed during the running of
// the interior point algorithm.  Each structure is allocated when an
// Dynamical Solver is initialized, and reused in each iteration.
//
class Dynamical_Solver
{
public:
  size_t total_iteration;
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
    duality_gap,                 // normalized difference of objectives
    external_step_size;          // the size of the step to be taken in the external parameters' space 


  El::Matrix<El::BigFloat> grad_BFGS, hess_BFGS, hess_BFGS_pp, grad_mixed;

  El::Matrix<El::BigFloat> hess_Exact;

  El::Matrix<El::BigFloat> grad_withlog;
  El::Matrix<El::BigFloat> grad_withoutlog;

  El::BigFloat lag_shifted;
  bool findMinimumQ;

  bool hess_BFGS_updateQ;

  El::BigFloat p_step, d_step;

  El::BigFloat prev_p_step, prev_d_step;

  El::BigFloat mulogdetX;

  std::vector<std::pair<El::BigFloat, El::BigFloat>> beta_scan_recorder;

  El::BigFloat final_beta;

  El::Matrix<El::BigFloat> specified_ext_param;
  bool specified_ext_param_Q;



  //int intext_mode;
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
  El::BigFloat R_error; // maxAbs(R)

  int64_t current_generation;
  boost::optional<int64_t> backup_generation;

  Dynamical_Solver(const Dynamical_Solver_Parameters &dynamical_parameters,
             const Verbosity &verbosity,
             const bool &require_initial_checkpoint,
             const Block_Info &block_info, const El::Grid &grid,
             const size_t &dual_objective_b_height);

  Dynamical_Solver_Terminate_Reason
  run_dynamical(const Dynamical_Solver_Parameters &dynamical_parameters,
      const Verbosity &verbosity,
      const SDP &sdp, 
      const boost::property_tree::ptree &parameter_properties,
      const Block_Info &block_info, const El::Grid &grid,
      Timers &timers, bool &update_sdp,  El::Matrix<El::BigFloat> &extParamStep);


  void external_corrector_step(const Dynamical_Solver_Parameters &dynamical_parameters,
	  const SDP &old_sdp, const SDP &new_sdp,
	  const Block_Info &block_info, const El::Grid &grid, const std::size_t &total_psd_rows,
	  Block_Diagonal_Matrix &X_cholesky, Block_Diagonal_Matrix &Y_cholesky,
	  std::array<std::vector<std::vector<std::vector<El::DistMatrix<El::BigFloat>>>>, 2> &A_X_inv,
	  std::array<std::vector<std::vector<std::vector<El::DistMatrix<El::BigFloat>>>>, 2> &A_Y,
	  Block_Vector & dx, Block_Vector & dy,
	  Block_Diagonal_Matrix & dX, Block_Diagonal_Matrix & dY,
	  const Block_Diagonal_Matrix &schur_complement_cholesky,
	  const Block_Matrix &schur_off_diagonal,
	  const El::DistMatrix<El::BigFloat> &Q,
	  El::BigFloat &current_error_R, El::BigFloat &current_mu,
	  El::BigFloat &primal_step_length, El::BigFloat &dual_step_length, El::BigFloat &step_length_reduction,
	  Timers &timers);

  void compute_corrector_residue_shift(const Block_Info &block_info,
	  Block_Vector &primal_residue_p_0, Block_Vector &dual_residues_0, Block_Diagonal_Matrix &R_0,
	  Block_Vector &primal_residue_p, Block_Vector &dual_residues, Block_Diagonal_Matrix &R,
	  Block_Vector & dx, Block_Vector & dy, Block_Diagonal_Matrix & dX, Block_Diagonal_Matrix & dY,
	  const SDP &d_sdp);

  void external_corrector_run(const Dynamical_Solver_Parameters &dynamical_parameters,
	  const SDP &new_sdp,
	  const Block_Info &block_info, const El::Grid &grid,
	  Block_Vector & dx, Block_Vector & dy, Block_Diagonal_Matrix & dX, Block_Diagonal_Matrix & dY,
	  El::BigFloat &primal_step_length, El::BigFloat &dual_step_length, El::BigFloat &step_length_reduction,
	  Timers &timers);

  void dynamical_step(
    const Dynamical_Solver_Parameters &dynamical_parameters,      
    const std::size_t &total_psd_rows,
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
    El::BigFloat &beta, El::BigFloat &primal_step_length, 
    El::BigFloat &dual_step_length, bool &terminate_now, Timers &timers,
    bool &update_sdp, bool &find_zeros, El::Matrix<El::BigFloat> &external_step); 

  void save_solver_state(const Dynamical_Solver_Parameters &dynamical_parameters,
	  const Block_Info &block_info,
	  Block_Vector & dx, Block_Vector & dy,
	  Block_Diagonal_Matrix & dX, Block_Diagonal_Matrix & dY,
	  const Block_Diagonal_Matrix &schur_complement_cholesky,
	  const Block_Matrix &schur_off_diagonal,
	  const El::DistMatrix<El::BigFloat> &Q);

  void save_checkpoint(
	  Block_Vector & dx, Block_Vector & dy,
	  Block_Diagonal_Matrix & dX, Block_Diagonal_Matrix & dY,
	  const boost::filesystem::path &checkpoint_directory,
	  const Verbosity &verbosity,
	  const boost::property_tree::ptree &parameter_properties);

  void
	  save_checkpoint(const boost::filesystem::path &checkpoint_directory,
		  const Verbosity &verbosity,
		  const boost::property_tree::ptree &parameter_properties);
  bool
	  load_checkpoint(const boost::filesystem::path &checkpoint_directory,
		  const Block_Info &block_info, const Verbosity &verbosity,
		  const bool &require_initial_checkpoint);


  void compute_external_dxdydXdY(
	  const bool &is_primal_and_dual_feasible,
	  const Dynamical_Solver_Parameters &dynamical_parameters,
	  const Block_Info &block_info,
	  const SDP &sdp, const El::Grid &grid,
	  const Block_Diagonal_Matrix &X_cholesky,
	  const Block_Diagonal_Matrix &Y_cholesky,
	  Timers &timers,

	  Block_Vector & internal_dx, Block_Vector & internal_dy,
	  Block_Vector & dx, Block_Vector & dy,
	  Block_Diagonal_Matrix & dX, Block_Diagonal_Matrix & dY,
	  Block_Diagonal_Matrix & R,

	  El::Matrix<El::BigFloat> & external_step,
	  std::vector<std::pair<Block_Vector, Block_Vector>> & Delta_xy,
	  El::BigFloat &primal_step_length, El::BigFloat &dual_step_length,
	  El::BigFloat &step_length_reduction
	  );

  void internal_step(
	  const Dynamical_Solver_Parameters &dynamical_parameters,
	  const std::size_t &total_psd_rows,
	  const Block_Info &block_info,
	  const SDP &sdp, const El::Grid &grid,
	  const Block_Diagonal_Matrix &X_cholesky,
	  const Block_Diagonal_Matrix &Y_cholesky,
	  Timers &timers,

	  const Block_Diagonal_Matrix &schur_complement_cholesky,
	  const Block_Matrix &schur_off_diagonal,
	  const El::DistMatrix<El::BigFloat> &Q,

	  Block_Vector & dx, Block_Vector & dy,
	  Block_Diagonal_Matrix & dX, Block_Diagonal_Matrix & dY,
	  Block_Diagonal_Matrix & R,
	  Block_Vector &grad_x, Block_Vector &grad_y,

	  const Block_Vector &primal_residue_p, El::BigFloat &mu,

	  const bool &is_primal_and_dual_feasible,
	  El::BigFloat &beta,
	  El::BigFloat &primal_step_length, El::BigFloat &dual_step_length,
	  El::BigFloat &step_length_reduction
  );

  void execute_step(
	  Block_Vector & dx, Block_Vector & dy,
	  Block_Diagonal_Matrix & dX, Block_Diagonal_Matrix & dY,
	  El::BigFloat &primal_step_length,
	  El::BigFloat &dual_step_length
  );


  void compute_dXdY(
	  const bool &is_primal_and_dual_feasible,
	  const Block_Info &block_info,
	  const SDP &sdp, const El::Grid &grid,
	  const Block_Diagonal_Matrix &X_cholesky,
	  const Block_Diagonal_Matrix &Y_cholesky,
	  Timers &timers,

	  Block_Vector & dx, Block_Vector & dy,
	  Block_Diagonal_Matrix & dX, Block_Diagonal_Matrix & dY,
	  Block_Diagonal_Matrix & R,
	  El::BigFloat &primal_step_length, El::BigFloat &dual_step_length,
	  El::BigFloat &step_length_reduction
  );

  void compute_dXdY(
	  const bool &is_primal_and_dual_feasible,
	  const Block_Info &block_info,
	  const SDP &sdp, const El::Grid &grid,
	  const Block_Diagonal_Matrix &X_cholesky,
	  const Block_Diagonal_Matrix &Y_cholesky,
	  Timers &timers,

	  Block_Vector & dx, Block_Vector & dy,
	  Block_Diagonal_Matrix & dX, Block_Diagonal_Matrix & dY,
	  Block_Diagonal_Matrix & R
  );

  void compute_external_dxdy(
	  const Dynamical_Solver_Parameters &dynamical_parameters,

	  Block_Vector & internal_dx, Block_Vector & internal_dy,
	  Block_Vector & dx, Block_Vector & dy,
	  Block_Diagonal_Matrix & R,

	  El::Matrix<El::BigFloat> & external_step,
	  std::vector<std::pair<Block_Vector, Block_Vector>> & Delta_xy
  );

  void strategy_hess_BFGS(const Dynamical_Solver_Parameters &dynamical_parameters, int n_external_parameters,
	  bool lowest_mu_Q,
	  El::Matrix<El::BigFloat> & grad_p, El::Matrix<El::BigFloat> & grad_mixed, El::Matrix<El::BigFloat> & grad_corrected,
	  El::Matrix<El::BigFloat> & Lpu, El::BigFloat & mu,
	  El::Matrix<El::BigFloat> & hess_pp, El::Matrix<El::BigFloat> & hess_mixed, El::Matrix<El::BigFloat> & hess_Exact,

	  El::Matrix<El::BigFloat> & prev_BFGS, El::Matrix<El::BigFloat> & prev_BFGS_pp, El::Matrix<El::BigFloat> & prev_step, El::Matrix<El::BigFloat> & prev_grad,
	  El::Matrix<El::BigFloat> & hess_BFGS_lowest_mu
  );


  void strategy_findboundary_extstep(
	  const Block_Info &block_info,
	  const Block_Diagonal_Matrix &X_cholesky,
	  const std::size_t &total_psd_rows,
	  const int dim_y,
	  const El::BigFloat & mu,
	  const El::BigFloat & beta,

	  int n_external_parameters,
	  const Dynamical_Solver_Parameters &dynamical_parameters,
	  bool lowest_mu_Q,

	  El::Matrix<El::BigFloat> & grad_corrected,
	  El::Matrix<El::BigFloat> & grad_p,
	  El::Matrix<El::BigFloat> & grad_mixed,

	  El::Matrix<El::BigFloat> & external_step
  );

  void strategy_update_grad_BFGS(
	  El::BigFloat &primal_step_length, El::BigFloat &dual_step_length,
	  El::Matrix<El::BigFloat> & grad_p, El::Matrix<El::BigFloat> & grad_mixed, El::Matrix<El::BigFloat> & grad_BFGS
  );

  void internal_step_corrector_iteration_centering(
	  const Dynamical_Solver_Parameters &dynamical_parameters,
	  const std::size_t &total_psd_rows,
	  const Block_Info &block_info,
	  const SDP &sdp, const El::Grid &grid,
	  const Block_Diagonal_Matrix &X_cholesky,
	  const Block_Diagonal_Matrix &Y_cholesky,
	  Timers &timers,

	  const Block_Diagonal_Matrix &schur_complement_cholesky,
	  const Block_Matrix &schur_off_diagonal,
	  const El::DistMatrix<El::BigFloat> &Q,

	  Block_Vector & dx, Block_Vector & dy,
	  Block_Diagonal_Matrix & dX, Block_Diagonal_Matrix & dY,
	  Block_Diagonal_Matrix & R,
	  Block_Vector &grad_x, Block_Vector &grad_y,

	  const Block_Vector &primal_residue_p, El::BigFloat &mu,

	  const bool &is_primal_and_dual_feasible,
	  El::BigFloat &beta,
	  El::BigFloat &primal_step_length, El::BigFloat &dual_step_length,
	  El::BigFloat &step_length_reduction
  );

  El::BigFloat finite_mu_navigator(
	  const Block_Info &block_info,
	  const Block_Diagonal_Matrix &X_cholesky,
	  const std::size_t &total_psd_rows,
	  const int dim_y,
	  const El::BigFloat & mu,
	  const El::BigFloat & beta,
	  const Dynamical_Solver_Parameters &dynamical_parameters);

  void dynamical_step_scan_beta(const Dynamical_Solver_Parameters &dynamical_parameters,
	  const El::Matrix<El::BigFloat> & eplus, const El::Matrix<El::BigFloat> & eminus, const El::Matrix<El::BigFloat> & esum,
	  const std::vector<std::pair<Block_Vector, Block_Vector>> & H_xp, std::vector<std::pair<Block_Vector, Block_Vector>> & Delta_xy,
	  Block_Vector & internal_dx, Block_Vector & internal_dy,
	  const El::Matrix<El::BigFloat> & grad_withoutlog,
	  const El::Matrix<El::BigFloat> & grad_withlog,
	  El::Matrix<El::BigFloat> & grad_p, El::Matrix<El::BigFloat> & grad_mixed, El::Matrix<El::BigFloat> & grad_corrected,
	  El::Matrix<El::BigFloat> & Lpu, El::BigFloat & mu,
	  El::Matrix<El::BigFloat> & hess_pp, El::Matrix<El::BigFloat> & hess_mixed, El::Matrix<El::BigFloat> & hess_Exact,

	  const std::size_t &total_psd_rows,
	  const Block_Info &block_info,
	  const SDP &sdp, const El::Grid &grid,
	  const Block_Diagonal_Matrix &X_cholesky,
	  const Block_Diagonal_Matrix &Y_cholesky,
	  Timers &timers,

	  const Block_Diagonal_Matrix &schur_complement_cholesky,
	  const Block_Matrix &schur_off_diagonal,
	  const El::DistMatrix<El::BigFloat> &Q,

	  Block_Vector & dx, Block_Vector & dy,
	  Block_Diagonal_Matrix & dX, Block_Diagonal_Matrix & dY,
	  Block_Diagonal_Matrix & R,
	  Block_Vector &grad_x, Block_Vector &grad_y,

	  const Block_Vector &primal_residue_p,

	  const bool &is_primal_and_dual_feasible,
	  El::BigFloat &beta,
	  El::BigFloat &primal_step_length, El::BigFloat &dual_step_length,
	  El::BigFloat &step_length_reduction,

	  int n_external_parameters,
	  El::Matrix<El::BigFloat> & prev_grad,
	  El::Matrix<El::BigFloat> & prev_step,
	  El::Matrix<El::BigFloat> & prev_BFGS,
	  El::Matrix<El::BigFloat> & external_step,
	  Block_Diagonal_Matrix & dX_best, Block_Diagonal_Matrix & dY_best,
	  Block_Vector & dx_best, Block_Vector & dy_best,
	  El::BigFloat & primal_step_length_best, El::BigFloat & dual_step_length_best, El::BigFloat & beta_best,
	  El::Matrix<El::BigFloat> & hess_BFGS_best, El::Matrix<El::BigFloat> & grad_mixed_best,
	  bool & update_sdp,

	  bool first_beta_Q
  );

};
