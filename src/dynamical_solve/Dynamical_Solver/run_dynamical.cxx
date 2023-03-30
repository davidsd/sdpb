#include "../Dynamical_Solver.hxx"
#include "../../sdp_solve.hxx"
#include "../../sdp_solve/SDP_Solver.hxx"
#include <iostream>
#include <boost/filesystem/fstream.hpp>

//Functions from SDP_Solver
void cholesky_decomposition(const Block_Diagonal_Matrix &A,
                            Block_Diagonal_Matrix &L);

void print_header_dynamical(const Verbosity &verbosity);
void print_header_dynamical_new(const Verbosity &verbosity);

void print_iteration(
  const int &iteration, const El::BigFloat &mu,
  const El::BigFloat &primal_step_length, const El::BigFloat &dual_step_length,
  const El::BigFloat &beta, const Dynamical_Solver &dynamical_solver,
  const std::chrono::time_point<std::chrono::high_resolution_clock>
    &solver_start_time,
  const Verbosity &verbosity);

void print_iteration_new(
	const int &iteration, const El::BigFloat &mu,
	const El::BigFloat &primal_step_length, const El::BigFloat &dual_step_length,
	const El::BigFloat &beta, const Dynamical_Solver &dynamical_solver,
	const std::chrono::time_point<std::chrono::high_resolution_clock>
	&solver_start_time,
	const Verbosity &verbosity);

void compute_objectives(const SDP &sdp, const Block_Vector &x,
                        const Block_Vector &y, El::BigFloat &primal_objective,
                        El::BigFloat &dual_objective,
                        El::BigFloat &duality_gap, Timers &timers);

void compute_bilinear_pairings(
  const Block_Info &block_info, const Block_Diagonal_Matrix &X_cholesky,
  const Block_Diagonal_Matrix &Y,
  const std::vector<El::DistMatrix<El::BigFloat>> &bases_blocks,
  std::array<std::vector<std::vector<std::vector<El::DistMatrix<El::BigFloat>>>>,
             2> &A_X_inv,
  std::array<std::vector<std::vector<std::vector<El::DistMatrix<El::BigFloat>>>>,
             2> &A_Y,
  Timers &timers);


void compute_dual_residues_and_error(
  const Block_Info &block_info, const SDP &sdp, const Block_Vector &y,
  const std::array<
    std::vector<std::vector<std::vector<El::DistMatrix<El::BigFloat>>>>, 2>
    &A_Y,
  Block_Vector &dual_residues, El::BigFloat &dual_error, Timers &timers);

void compute_primal_residues_and_error_P_Ax_X(
  const Block_Info &block_info, const SDP &sdp, const Block_Vector &x,
  const Block_Diagonal_Matrix &X, Block_Diagonal_Matrix &primal_residues,
  El::BigFloat &primal_error_P, Timers &timers);

void compute_primal_residues_and_error_p_b_Bx(const Block_Info &block_info,
                                              const SDP &sdp,
                                              const Block_Vector &x,
                                              Block_Vector &primal_residue_p,
                                              El::BigFloat &primal_error_p);


void compute_R_error(const std::size_t &total_psd_rows, const Block_Diagonal_Matrix &X, const Block_Diagonal_Matrix &Y, El::BigFloat & R_error, Timers &timers);

extern int max_climbing;

// subroutines to decide weather to update sdps in run_dynamical, before entering dynamical_step
void compute_update_sdp(
  const Dynamical_Solver_Parameters &parameters, const El::BigFloat &primal_error,
  const El::BigFloat &dual_error, const El::BigFloat &duality_gap,
  const El::BigFloat &primal_step_length, const El::BigFloat &dual_step_length,
  const int &iteration,
  const std::chrono::time_point<std::chrono::high_resolution_clock>
    &solver_start_time,
  bool &is_primal_and_dual_feasible,
  bool &update_sdp);

//Same as the function in SDP_Solver, adapted for the Dynamical_Solver_Terminate_Reason class
void compute_feasible_and_termination(
  const Solver_Parameters &parameters, const El::BigFloat &primal_error,
  const El::BigFloat &dual_error, const El::BigFloat &duality_gap,
  const El::BigFloat &primal_step_length, const El::BigFloat &dual_step_length,
  const int &iteration,
  const std::chrono::time_point<std::chrono::high_resolution_clock>
    &solver_start_time,
  bool &is_primal_and_dual_feasible,
  Dynamical_Solver_Terminate_Reason &terminate_reason, bool &terminate_now);

// checkpoint_save will use this variable. It's very bad we pass this variable as argument everywhere.
// I will just define a global variable. This is temporary fix. We should make it member variable.
boost::property_tree::ptree parameter_properties_save;
Verbosity verbosity_save;

Dynamical_Solver_Terminate_Reason
Dynamical_Solver::run_dynamical(const Dynamical_Solver_Parameters &dynamical_parameters,
                const Verbosity &verbosity,
                const SDP &sdp, 
                const boost::property_tree::ptree &parameter_properties,
                const Block_Info &block_info, 
                const El::Grid &grid, Timers &timers, 
                bool &update_sdp, El::Matrix<El::BigFloat> &extParamStep)
{
	parameter_properties_save = parameter_properties;
	verbosity_save = verbosity;
	max_climbing = dynamical_parameters.max_climbing;

  Dynamical_Solver_Terminate_Reason terminate_reason(
    Dynamical_Solver_Terminate_Reason::MaxIterationsExceeded);
  auto &solver_timer(timers.add_and_start("Solver runtime"));
  auto &initialize_timer(timers.add_and_start("run.initialize"));

  El::BigFloat primal_step_length(0), dual_step_length(0);

  Block_Diagonal_Matrix X_cholesky(X), Y_cholesky(X);

  std::array<
    std::vector<std::vector<std::vector<El::DistMatrix<El::BigFloat>>>>, 2>
    A_X_inv;

  std::array<
    std::vector<std::vector<std::vector<El::DistMatrix<El::BigFloat>>>>, 2>
    A_Y;

  if (dynamical_parameters.printMore)
	  print_header_dynamical_new(verbosity);
  else
	  print_header_dynamical(verbosity);

  auto psd_sizes(block_info.psd_matrix_block_sizes());
  std::size_t total_psd_rows(
    std::accumulate(psd_sizes.begin(), psd_sizes.end(), size_t(0)));

  initialize_timer.stop();
  auto last_checkpoint_time(std::chrono::high_resolution_clock::now());
  size_t iteration  = 1; 
  bool find_zeros = false;
  El::BigFloat mu;
  for(; ; ++iteration)
    {
      El::byte checkpoint_now(
        std::chrono::duration_cast<std::chrono::seconds>(
          std::chrono::high_resolution_clock::now() - last_checkpoint_time)
          .count()
        >= dynamical_parameters.solver_parameters.checkpoint_interval);
      // Time varies between cores, so follow the decision of the root.
      El::mpi::Broadcast(checkpoint_now, 0, El::mpi::COMM_WORLD);
      if(checkpoint_now == true)
        {
          save_checkpoint(dynamical_parameters.solver_parameters.checkpoint_out, verbosity,
                          parameter_properties);
          last_checkpoint_time = std::chrono::high_resolution_clock::now();
        }
      compute_objectives(sdp, x, y, primal_objective, dual_objective,
                         duality_gap, timers);

      auto &cholesky_decomposition_timer(
        timers.add_and_start("run.choleskyDecomposition"));
      cholesky_decomposition(X, X_cholesky);
      cholesky_decomposition(Y, Y_cholesky);
      cholesky_decomposition_timer.stop();

      compute_bilinear_pairings(block_info, X_cholesky, Y, sdp.bases_blocks,
                                A_X_inv, A_Y, timers);

      compute_dual_residues_and_error(block_info, sdp, y, A_Y, dual_residues,
                                      dual_error, timers);
      compute_primal_residues_and_error_P_Ax_X(
        block_info, sdp, x, X, primal_residues, primal_error_P, timers);

      // use y to set the sizes of primal_residue_p.  The data is
      // overwritten in compute_primal_residues_and_error_p.
      Block_Vector primal_residue_p(y);
      compute_primal_residues_and_error_p_b_Bx(
        block_info, sdp, x, primal_residue_p, primal_error_p);

	  compute_R_error(total_psd_rows, X, Y, R_error, timers);

      bool terminate_now, is_primal_and_dual_feasible;
      compute_feasible_and_termination(
        dynamical_parameters.solver_parameters, primal_error(), dual_error, duality_gap,
        primal_step_length, dual_step_length, iteration,
        solver_timer.start_time, is_primal_and_dual_feasible, terminate_reason,
        terminate_now);
      if(terminate_now)
        {
          break;
        }

      //El::BigFloat mu; 
      El::BigFloat beta_predictor;
       {
         external_step_size = 0; 
         dynamical_step(dynamical_parameters, total_psd_rows, is_primal_and_dual_feasible, block_info,
                      sdp, grid, 
                      X_cholesky, Y_cholesky, A_X_inv, A_Y,  primal_residue_p, 
                      mu, beta_predictor, primal_step_length, dual_step_length, 
                      terminate_now, timers, update_sdp, find_zeros, extParamStep);
		 final_beta = beta_predictor;
       } 
      if (terminate_now)
        {
          terminate_reason
            = Dynamical_Solver_Terminate_Reason::MaxComplementarityExceeded;
          break;
        }


	  if (dynamical_parameters.printMore)
		  print_iteration_new(iteration, mu, primal_step_length, dual_step_length,
			  beta_predictor, *this, solver_timer.start_time,
			  verbosity);
	  else
		  print_iteration(iteration, mu, primal_step_length, dual_step_length,
			  beta_predictor, *this, solver_timer.start_time,
			  verbosity);


     if (update_sdp)
       {
         terminate_reason
            = Dynamical_Solver_Terminate_Reason::UpdateSDPs;
         break;
       }
    
    }

  if (El::mpi::Rank() == 0) std::cout << "hess_BFGS_updateQ = " << hess_BFGS_updateQ << "\n" << std::flush;

  if (terminate_reason != Dynamical_Solver_Terminate_Reason::UpdateSDPs)
     { El::Zero(extParamStep);}
  total_iteration = dynamical_parameters.total_iterations + iteration; 
  solver_timer.stop();
  return terminate_reason;
}
