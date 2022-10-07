#include <unistd.h>

#include <chrono>

#include <vector>
#include <boost/filesystem.hpp>
#include <boost/algorithm/string.hpp>
#include <iostream>
#include "../../../sdp_solve.hxx"
#include "../../../dynamical_solve.hxx"
#include "../../../sdp_read.hxx"
#include "imports.hxx"
#include "../../../approx_objective/Approx_Objective.hxx"

//Compute dx and dy of the central sdp as the standard sdp_solver does in predictor phase. 
//Correspond to - H^-1_xx Del_p L_mu in Eq(13).
//Return: void.
//Update dx, dy. 
void internal_predictor_direction(
	const Block_Info &block_info, const SDP &sdp, const Dynamical_Solver &solver,
	const Block_Diagonal_Matrix &schur_complement_cholesky,
	const Block_Matrix &schur_off_diagonal,
	const Block_Diagonal_Matrix &X_cholesky, const El::BigFloat beta,
	const El::BigFloat &mu, const Block_Vector &primal_residue_p,
	const El::DistMatrix<El::BigFloat> &Q, Block_Vector &grad_x, Block_Vector &grad_y, Block_Vector &dx, Block_Vector &dy, Block_Diagonal_Matrix &R);

void internal_corrector_direction(
	const Block_Info &block_info, const SDP &sdp, const Dynamical_Solver &solver,
	const Block_Diagonal_Matrix &schur_complement_cholesky,
	const Block_Matrix &schur_off_diagonal,
	const Block_Diagonal_Matrix &X_cholesky, const El::BigFloat beta,
	const El::BigFloat &mu, const Block_Vector &primal_residue_p,
	const El::DistMatrix<El::BigFloat> &Q, Block_Vector &grad_x, Block_Vector &grad_y, Block_Vector &dx, Block_Vector &dy, Block_Diagonal_Matrix &dX, Block_Diagonal_Matrix &dY, Block_Diagonal_Matrix &R);

void internal_predictor_direction_dxdydXdY(
	const Block_Info &block_info, const SDP &sdp, const Dynamical_Solver &solver,
	const Block_Diagonal_Matrix &schur_complement_cholesky,
	const Block_Matrix &schur_off_diagonal,
	const Block_Diagonal_Matrix &X_cholesky, const El::BigFloat beta,
	const El::BigFloat &mu, const Block_Vector &primal_residue_p,
	const El::DistMatrix<El::BigFloat> &Q, Block_Vector &grad_x, Block_Vector &grad_y, Block_Vector &dx, Block_Vector &dy, Block_Diagonal_Matrix &dX, Block_Diagonal_Matrix &dY, Block_Diagonal_Matrix &R);


//Compute the linear difference of the lagrangians of two sdps 
//The same as the calculation done by Approx_Objective
El::BigFloat compute_delta_lag(const SDP &d_sdp, const Dynamical_Solver &solver);

El::BigFloat compute_lag(const El::BigFloat mu, const Block_Diagonal_Matrix &X_cholesky, const Dynamical_Solver &solver);

//Given delta_p(sdp) , compute the (delta_x, delta_y) = H_xx^-1 H_xp delta_p
//as shown in Eq(15). 
//Return: void. 
//Update: hess_xp = H_xp = (RHS(p1)/alpha, RHS(p2)/alpha, ... ), stored to compute the second term on the LHS of Eq(13) 
//        delta_x_y = - H^-1_xx H_xp = - H^-1_xx hess_xp 
void mixed_hess(
	const Block_Info &block_info, const SDP &d_sdp,
	const Block_Vector &x, const Block_Vector &y,
	const Block_Diagonal_Matrix &schur_complement_cholesky,
	const Block_Matrix &schur_off_diagonal,
	const El::DistMatrix<El::BigFloat> &Q,
	const El::BigFloat &alpha,
	std::vector<std::pair<Block_Vector, Block_Vector>> &hess_xp,
	std::vector<std::pair<Block_Vector, Block_Vector>> &delta_x_y);

//Given objectives on a grid in the external parameter (p) space,
//where ePlus correspond to (+e_i), eMinus corresponds to (-e_i) and eSum corresponds to (e_i+e_j) directions respectively, 
//Compute H_pp and Del_p(L_mu) to solve Eq(13).
//Return: void.
//Update: hess = H_pp ,
//        grad = Del_p(L_mu).
void external_grad_hessian(const El::Matrix<El::BigFloat> &ePlus,
	const El::Matrix<El::BigFloat> &eMinus,
	const El::Matrix<El::BigFloat> &eSum,
	const El::BigFloat &alpha,
	El::Matrix<El::BigFloat> &grad,
	El::Matrix<El::BigFloat> &hess);


void external_window(const std::vector<El::BigFloat> &bmax,
	const std::vector<El::BigFloat> &bmin,
	const std::vector<El::BigFloat> &external_coord,
	const int n_parameters,
	El::Matrix<El::BigFloat> &step)
{
	El::Print(step);
	for (int i = 0; i < n_parameters; i++) {
		step(i) = El::Min(El::Max(bmin.at(i) - external_coord.at(i), step(i)), bmax.at(i) - external_coord.at(i));
	}
	El::Print(step);
}

//void external_grad(const El::Matrix<El::BigFloat> &ePlus,
//	const El::Matrix<El::BigFloat> &eMinus,
//	const El::BigFloat &alpha,
//	El::Matrix<El::BigFloat> &grad);

//This function is the same as in the usual SDP_Solver, 
//modified to take the argument 'solver' of the right type. 
//Currently not used in the dynamical_step 
void compute_search_direction(
	const Block_Info &block_info, const SDP &sdp, const Dynamical_Solver &solver,
	const Block_Diagonal_Matrix &schur_complement_cholesky,
	const Block_Matrix &schur_off_diagonal,
	const Block_Diagonal_Matrix &X_cholesky, const El::BigFloat beta,
	const El::BigFloat &mu, const Block_Vector &primal_residue_p,
	const bool &is_corrector_phase, const El::DistMatrix<El::BigFloat> &Q,
	Block_Vector &dx, Block_Diagonal_Matrix &dX, Block_Vector &dy,
	Block_Diagonal_Matrix &dY);

void compute_find_zeros(
	const El::BigFloat &duality_gap, const El::BigFloat &primal_objective,
	bool &find_zeros);


void find_zero_step(const El::BigFloat &thresh, const int &max_it, const El::BigFloat &step_size_max,
	El::Matrix<El::BigFloat> &hess, const El::Matrix<El::BigFloat> &grad, bool &find_zeros,
	El::Matrix<El::BigFloat> &external_step, const El::BigFloat &quadratic_const);


// Centering parameter \beta_p for the predictor step
El::BigFloat predictor_centering_parameter_V2(const Solver_Parameters &parameters,
	const bool is_primal_dual_feasible)
{
	return is_primal_dual_feasible ? parameters.feasible_centering_parameter
		: parameters.infeasible_centering_parameter;
}


// << std::setprecision(100) 
void print_matrix(const El::Matrix<El::BigFloat> & matrix)
{
	int dim_width = matrix.Width();
	int dim_height = matrix.Height();

	for (int i = 0; i < dim_height; i++)
	{
		for (int j = 0; j < dim_width; j++)if (El::mpi::Rank() == 0) std::cout << matrix(i, j) << " ";
		if (El::mpi::Rank() == 0) std::cout << "\n";
	}
	return;
}

void print_matrix(const El::Matrix<El::BigFloat> & matrix, int dim)
{
	int dim_width = std::min(matrix.Width(), dim);
	int dim_height = std::min(matrix.Height(), dim);

	for (int i = 0; i < dim_height; i++)
	{
		for (int j = 0; j < dim_width; j++)if (El::mpi::Rank() == 0) std::cout << matrix(i, j) << " ";
		if (El::mpi::Rank() == 0) std::cout << "\n";
	}
	return;
}

void print_vector(const El::Matrix<El::BigFloat> & vec)
{
	for (int i = 0; i < vec.Height(); i++) if (El::mpi::Rank() == 0)std::cout << vec(i, 0) << " ";
}


bool positivitize_matrix(El::Matrix<El::BigFloat> & matrix);

void BFGS_update_hessian(const int n_parameters,
	const El::Matrix<El::BigFloat> &grad_p_diff,
	const El::Matrix<El::BigFloat> &last_it_step,
	El::Matrix<El::BigFloat> &hess_bfgs
);

void read_sdp_grid(
	const Dynamical_Solver_Parameters &dynamical_parameters,
	const Block_Info &block_info,
	const SDP &sdp, const El::Grid &grid,
	Timers &timers,

	Block_Diagonal_Matrix & schur_complement_cholesky,
	Block_Matrix & schur_off_diagonal,

	El::DistMatrix<El::BigFloat> & Q,

	Block_Vector & x,
	Block_Vector & y,

	const Block_Diagonal_Matrix &X_cholesky,

	int n_external_parameters,

	El::Matrix<El::BigFloat> & eplus,
	El::Matrix<El::BigFloat> & eminus,
	El::Matrix<El::BigFloat> & esum,

	El::Matrix<El::BigFloat> & Lpu,

	std::vector<std::pair<Block_Vector, Block_Vector>> & H_xp,
	std::vector<std::pair<Block_Vector, Block_Vector>> & Delta_xy);


void compute_ellipse_boundary(const El::Matrix<El::BigFloat> &H,
	const El::Matrix<El::BigFloat> &g,
	El::BigFloat &f0,
	const El::Matrix<El::BigFloat> &search_direction,
	El::Matrix<El::BigFloat> & result_plus,
	El::Matrix<El::BigFloat> & result_minus,
	El::BigFloat & lambda);


bool BFGS_update_hessian(
	const El::Matrix<El::BigFloat> &grad_p_diff,
	const El::Matrix<El::BigFloat> &last_it_step,
	El::Matrix<El::BigFloat> &hess_bfgs,
	bool update_only_when_positive);

bool BFGS_partial_update_hessian(const El::BigFloat & reduction_factor,
	const El::Matrix<El::BigFloat> &y,
	const El::Matrix<El::BigFloat> &s,
	El::Matrix<El::BigFloat> &B);


void read_prev_grad_step_hess(const Dynamical_Solver_Parameters &dynamical_parameters,
	El::Matrix<El::BigFloat> & prev_grad, El::Matrix<El::BigFloat> & prev_step, El::Matrix<El::BigFloat> & prev_BFGS);

void compute_grad_p_grad_mixed_Hpp_Hmixed(const Dynamical_Solver_Parameters &dynamical_parameters,
	const El::Matrix<El::BigFloat> & eplus, const El::Matrix<El::BigFloat> & eminus, const El::Matrix<El::BigFloat> & esum,
	const std::vector<std::pair<Block_Vector, Block_Vector>> & H_xp, const std::vector<std::pair<Block_Vector, Block_Vector>> & Delta_xy,
	const Block_Vector & internal_dx, const Block_Vector & internal_dy,
	El::Matrix<El::BigFloat> & grad_p, El::Matrix<El::BigFloat> & grad_mixed,
	El::Matrix<El::BigFloat> & hess_pp, El::Matrix<El::BigFloat> & hess_mixed);


void save_load_beta_scan_best_state(Block_Vector & dx, Block_Vector & dy,
	Block_Diagonal_Matrix & dX, Block_Diagonal_Matrix & dY,
	El::BigFloat & primal_step_length, El::BigFloat & dual_step_length,
	El::BigFloat & beta,
	El::Matrix<El::BigFloat> & hess_BFGS, El::Matrix<El::BigFloat> & grad_mixed,

	Block_Vector & dx_best, Block_Vector & dy_best,
	Block_Diagonal_Matrix & dX_best, Block_Diagonal_Matrix & dY_best,
	El::BigFloat & primal_step_length_best, El::BigFloat & dual_step_length_best,
	El::BigFloat & beta_best,
	El::Matrix<El::BigFloat> & hess_BFGS_best, El::Matrix<El::BigFloat> & grad_mixed_best
);


void extrapolate_gradient_as_target_mu(
	const El::BigFloat & primal_step_length, const El::BigFloat & dual_step_length,
	const El::Matrix<El::BigFloat> & grad_p, const El::Matrix<El::BigFloat> & grad_mixed_best, El::Matrix<El::BigFloat> & grad_BFGS);



// A subroutine called by run_dynamical
// The external parameter step is passed by the argument 'external_step'
// Compute external_step using Eq (13)
// Compute (dx, dy, dX, dY) using Eq(12)
// Scale both external_step and (dx, dy, dX, dY) by the step length 


El::Matrix<El::BigFloat> external_step_save;
bool external_step_specified_Q = false;

El::Matrix<El::BigFloat> hess_BFGS_lowest_mu;

bool lowest_mu_Q = true;  // true : the solver hasn't been lifted along the local central path

/*  // this is parameter for the "best" Ising run

bool compute_ext_step_only_once = false;

int max_climbing = 1000;

bool update_hess_only_positive = false;

bool use_Lpu_mu_correction = false;

double PDsteplen_threshold = 0.3;

double navigator_shift = 0.2;

double dgap_near_optimal_threshold = 1e-20;

double beta_scan_start = -1;

*/


/*
bool compute_ext_step_only_once = true;

int max_climbing = 1000;

bool update_hess_only_positive = true;

bool use_Lpu_mu_correction = false;

double PDsteplen_threshold = 0.3;

double navigator_shift = 0.2;

double dgap_near_optimal_threshold = 1e-20;

double beta_scan_start = 0;  // if negative, scan starts with the value from predictor_centering_parameter_V3

*/



bool compute_ext_step_only_once = true;
bool recompute_ext_during_re_aiming = true;  // if true, ext-step will be recomputed during re-aiming at the lowest mu, even if compute_ext_step_only_once = true

int max_climbing = 2;

bool update_hess_only_positive = true;

bool use_Lpu_mu_correction = false;

double step_min_threshold = 0.1;
double step_max_threshold = 0.6;

double beta_scan_step = 0.1;

double navigator_shift = 0.2;

double dgap_near_optimal_threshold = 1e-20;

double beta_scan_start = -1;

bool use_gradp_for_BFGS_update = true;

double rescale_initial_hess = 1;

double BFGS_partial_update_reduction = -1; // if positive, hessian will be updated always positive based on the partial update logic.

bool compare_BFGS_gradient_at_same_mu = false;

double beta_for_mu_logdetX = 0; // if negative, it will use the beta from beta scan


El::BigFloat predictor_centering_parameter_V3(const Solver_Parameters &parameters,
	const El::BigFloat & mu, const El::BigFloat & mu_near_optimal_threshold)
{
	return mu < mu_near_optimal_threshold ? parameters.feasible_centering_parameter
		: parameters.infeasible_centering_parameter;
}

void Dynamical_Solver::dynamical_step(
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
	bool &update_sdp, bool &find_zeros, El::Matrix<El::BigFloat> &external_step)
{
	auto &step_timer(timers.add_and_start("run.step"));

	bool centering_mode(false);

	El::BigFloat delta_lambda(0);
	Block_Diagonal_Matrix schur_complement_cholesky(
		block_info.schur_block_sizes(), block_info.block_indices,
		block_info.num_points.size(), grid);

	Block_Matrix schur_off_diagonal(sdp.free_var_matrix);

	El::DistMatrix<El::BigFloat> Q(sdp.dual_objective_b.Height(),
		sdp.dual_objective_b.Height());

	initialize_schur_complement_solver(block_info, sdp, A_X_inv, A_Y, grid,
		schur_complement_cholesky,
		schur_off_diagonal, Q, timers);

	auto &frobenius_timer(
		timers.add_and_start("run.step.frobenius_product_symmetric"));
	mu = frobenius_product_symmetric(X, Y) / total_psd_rows;
	frobenius_timer.stop();
	if (mu > dynamical_parameters.solver_parameters.max_complementarity)
	{
		terminate_now = true;
		return;
	}

	//if (El::mpi::Rank() == 0) std::cout << "test A \n" << std::flush;

	// -----
	//Compute the predictor solution for (dx, dX, dy, dY)
	// -----
	auto &predictor_timer(
		timers.add_and_start("run.step.computeSearchDirection(betaPredictor)"));
	beta = predictor_centering_parameter_V3(dynamical_parameters.solver_parameters, mu*total_psd_rows, dgap_near_optimal_threshold);
	if (dynamical_parameters.total_iterations > 0 && beta_scan_start >= 0)
		beta = beta_scan_start;

	El::BigFloat step_length_reduction = dynamical_parameters.solver_parameters.step_length_reduction;
	if (mu*total_psd_rows < dgap_near_optimal_threshold)
	{
		//step_length_reduction = 0.9;
		beta = 0.1;
	}

	//std::cout << "rank=" << El::mpi::Rank() << " pid=" << getpid() << " check R_error=" << R_error << "\n" << std::flush;

	if (dynamical_parameters.updateSDP_dualityGapThreshold <= 0 &&
		dynamical_parameters.centeringRThreshold > 0 && R_error > dynamical_parameters.centeringRThreshold)
	{
		if (El::mpi::Rank() == 0)std::cout << "run centering steps with beta=1 \n";
		beta = 1;
		update_sdp = false;
		centering_mode = true;
	}

	//if (mu_direction_mode == 1)beta = 1;
	//if (mu_direction_mode == 2)beta = 2;

	//std::cout << "status : intext_mode=" << intext_mode << " beta=" << beta << " update_sdp=" << update_sdp <<  " \n";

	if (intext_mode == 1)
	{
		beta = 2;
		max_climbing--;
		lowest_mu_Q = false;
		intext_mode = 0;
	}

	// Internal_step: compute dx and dy for the central sdp as in compute_search_direction()      
	//                - H^-1_xx Del_x L_mu in Eq (12) and Eq(13)
	//                Notice that the negative sign has been accounted. 
	Block_Vector internal_dx(x), internal_dy(y);
	Block_Vector dx(internal_dx), dy(internal_dy);
	Block_Diagonal_Matrix dX(X), dY(Y);
	Block_Vector grad_x(internal_dx), grad_y(internal_dy);
	Block_Diagonal_Matrix R(X);

	Block_Diagonal_Matrix dX_best(X), dY_best(Y);
	Block_Vector dx_best(x), dy_best(y);
	El::BigFloat primal_step_length_best, dual_step_length_best, beta_best;
	El::Matrix<El::BigFloat> hess_BFGS_best, grad_mixed_best;

	//std::cout << "rank=" << El::mpi::Rank() << " pid=" << getpid() << " before for, update_sdp=" << update_sdp << "\n" << std::flush;

	int n_external_parameters = dynamical_parameters.n_external_parameters;
	std::vector<std::pair<Block_Vector, Block_Vector>> H_xp, Delta_xy;
	El::Matrix<El::BigFloat> eplus(n_external_parameters, 1), eminus(n_external_parameters, 1), esum(n_external_parameters, n_external_parameters);
	El::Matrix<El::BigFloat> Lpu(n_external_parameters, 1);

	El::Matrix<El::BigFloat> grad_p(n_external_parameters, 1); // Del_p L
	El::Matrix<El::BigFloat> hess_pp(n_external_parameters, n_external_parameters); //H_pp

	El::Matrix<El::BigFloat> hess_mixed(n_external_parameters, n_external_parameters); //H_px H^-1_xx H_xp in Eq(13).
	El::Matrix<El::BigFloat> grad_mixed(n_external_parameters, 1); //H_px (-internal_dx_dy) = H_px (H^-1_xx Del_x L_mu) in Eq (13).  


	if (update_sdp)
		read_sdp_grid(dynamical_parameters, block_info, sdp, grid, timers,
			schur_complement_cholesky, schur_off_diagonal, Q, x, y, X_cholesky, n_external_parameters,
			eplus, eminus, esum, Lpu,
			H_xp, Delta_xy);


	for (int stallrecovery_phase = 1; stallrecovery_phase <= 20; stallrecovery_phase++)
	{
		if (beta <= 10)
		{
			internal_predictor_direction_dxdydXdY(block_info, sdp, *this, schur_complement_cholesky,
				schur_off_diagonal, X_cholesky, beta,
				mu, primal_residue_p, Q, grad_x, grad_y, internal_dx, internal_dy, dX, dY, R);
			internal_corrector_direction(block_info, sdp, *this, schur_complement_cholesky,
				schur_off_diagonal, X_cholesky, beta,
				mu, primal_residue_p, Q, grad_x, grad_y, internal_dx, internal_dy, dX, dY, R);
		}
		else
		{
			internal_predictor_direction(block_info, sdp, *this, schur_complement_cholesky,
				schur_off_diagonal, X_cholesky, beta,
				mu, primal_residue_p, Q, grad_x, grad_y, internal_dx, internal_dy, R);
		}

		//if (El::mpi::Rank() == 0)std::cout << "test 3 \n" << std::flush;

		//Delta_xy = - H^-1_xx H_xp as the second term on the LHS of Eq(13)

		// Eq(13):
		// (H_pp - H_px H^-1_xx H_xp)      dp     =  - Del_p L   + H_px (H^-1_xx Del_x L)
		// (H_pp - (H_px H^-1_xx H_xp))    dp     =  - Del_p L   + H_px . (- internal_dx, - internal_dy)  
		// (hess_pp - hess_mixed)          dp     =  - grad_p    + grad_mixed 

		// External_step: compute hess_mixed := H_px H^-1_xx H_xp   , the second term on the LHS of Eq(13)
		//                        hess_pp    := H_pp                , the first term on the LHS of Eq(13)
		//                        grad_p     := Del_p L             , the first term on the RHS of Eq(13)
		//                        grad_mixed := H_px H^-1_xx Del_x L, the second term on the RHS of Eq(13)
		//                if criteria base on quantities, such as duality_gap, are met
		//                which is decided by function "compute_update_sdp()" called by run_dynamical()

		//if (El::mpi::Rank() == 0)std::cout << "test 4 \n" << std::flush;

		//std::cout << "rank=" << El::mpi::Rank() << " pid=" << getpid() << " update_sdp=" << update_sdp << " before if (update_sdp) \n" << std::flush;

		if (update_sdp)
		{
			compute_grad_p_grad_mixed_Hpp_Hmixed(dynamical_parameters, 
				eplus, eminus, esum, H_xp, Delta_xy, internal_dx, internal_dy,
				grad_p, grad_mixed, hess_pp, hess_mixed);

			// Eq(13):
			// (H_pp - H_px H^-1_xx H_xp)      dp     =  - Del_p L   + H_px (H^-1_xx Del_x L)
			// (H_pp - (H_px H^-1_xx H_xp))    dp     =  - Del_p L   + H_px . (- internal_dx, - internal_dy)  
			// (hess_pp - hess_mixed)          dp     =  - grad_p    + grad_mixed 

			// H_pp - H_px H^-1_xx H_xp, LHS of Eq(13)

			if (dynamical_parameters.use_exact_hessian)
			{
				hess_Exact = hess_pp;
				hess_Exact -= hess_mixed;
			}

			// - H_px (H^-1_xx Del_x L_mu) + Del_p L_mu, RHS of Eq(13) 

			El::Matrix<El::BigFloat> grad_corrected = grad_mixed;
			grad_corrected += grad_p;

			if (use_Lpu_mu_correction)
			{
				grad_corrected = Lpu;
				grad_corrected *= -mu;
				grad_corrected += grad_p;
			}

			/*
			if (use_gradp_for_BFGS_update)
				grad_BFGS = grad_p;
			else
				grad_BFGS = grad_corrected;
			*/

                        El::Zeros(hess_Exact, n_external_parameters, n_external_parameters);
			El::Matrix<El::BigFloat> prev_grad(n_external_parameters, 1), prev_step(n_external_parameters, 1), prev_BFGS(n_external_parameters, n_external_parameters);
			read_prev_grad_step_hess(dynamical_parameters, prev_grad, prev_step, prev_BFGS);

			// decide hess_BFGS
			if (dynamical_parameters.use_exact_hessian
				//|| (mu < 1e-8 && mu > 1e-9)
				)
			{
				if (positivitize_matrix(hess_Exact) == false)
					if (El::mpi::Rank() == 0) std::cout << "flip the sign of hessian\n" << std::flush;

				hess_BFGS = hess_Exact;
			}
			else 
			{
				if (El::mpi::Rank() == 0)std::cout << "using BFGS: " << '\n' << std::flush;

				if (dynamical_parameters.use_Hmixed_for_BFGS) // I will let simpleboot control whether hess_mixed will be used for hess_BFGS
				{
					El::Zeros(hess_BFGS, n_external_parameters, n_external_parameters);
					hess_BFGS -= hess_mixed;
					positivitize_matrix(hess_BFGS);
					hess_BFGS *= El::BigFloat(rescale_initial_hess);
				}
				else
				{
					hess_BFGS = prev_BFGS;
				}

				if (lowest_mu_Q == false) hess_BFGS = hess_BFGS_lowest_mu;

				// update hess_BFGS
				if (dynamical_parameters.total_iterations > 0 && lowest_mu_Q)
				{
					El::Matrix<El::BigFloat> grad_diff;
					if (use_gradp_for_BFGS_update)
						grad_diff = grad_p;
					else
						grad_diff = grad_corrected;

					grad_diff -= prev_grad;

					if (BFGS_partial_update_reduction > 0)
					{
						bool exact_update_Q = BFGS_partial_update_hessian(El::BigFloat(BFGS_partial_update_reduction),
							grad_diff, prev_step, hess_BFGS);
						if (!exact_update_Q)
							if (El::mpi::Rank() == 0)
								std::cout << "New hessian is non-positive. partial update. \n" << std::flush;
					}
					else
					{
						bool flippedQ = BFGS_update_hessian(grad_diff, prev_step, hess_BFGS, update_hess_only_positive);
						if (flippedQ)
							if (update_hess_only_positive)
							{
								if (El::mpi::Rank() == 0) std::cout << "New hessian is non-positive. Rollback to prev_hess. \n" << std::flush;
							}
							else
							{
								if (El::mpi::Rank() == 0) std::cout << "New hessian is non-positive (still updated). \n" << std::flush;
							}

						if (flippedQ == false || update_hess_only_positive == false)
							hess_BFGS_updateQ = true;
					}

					hess_BFGS_lowest_mu = hess_BFGS;

					//if (El::mpi::Rank() == 0) std::cout << "hess_BFGS_lowest_mu saved : " << std::flush;
					//print_matrix(hess_BFGS_lowest_mu);
				}
			}

			//if (El::mpi::Rank() == 0) std::cout << "hess_BFGS_lowest_mu : " << std::flush;
			//print_matrix(hess_BFGS_lowest_mu);

			//if (El::mpi::Rank() == 0) std::cout << "hess_BFGS : " << std::flush;
			//print_matrix(hess_BFGS);

			external_step = grad_corrected;
			external_step *= (-1);
			El::LinearSolve(hess_BFGS, external_step);
			external_step_size = El::Nrm2(external_step);

			if (dynamical_parameters.find_boundary && El::Abs(duality_gap) < 0.9 && dynamical_parameters.total_iterations > 0)
			{
				El::BigFloat lag0 = compute_lag(0, X_cholesky, *this);

				El::BigFloat beta_lag;
				if (beta_for_mu_logdetX < 0)
				{
					beta_lag = beta;
				}
				else
				{
					beta_lag = El::BigFloat(beta_for_mu_logdetX);
				}

				El::BigFloat lag = compute_lag(beta_lag*mu, X_cholesky, *this);

				El::BigFloat lag_shifted = lag;
				lag_shifted -= total_psd_rows * mu * El::BigFloat(navigator_shift);

				if (El::mpi::Rank() == 0)
				{
					std::cout << "lag_shifted  = " << lag_shifted << '\n' << std::flush;
					std::cout << "logdetX  = " << lag - lag0 << '\n' << std::flush;
				}

				//find_zero_step(10^(-10), 100, dynamical_parameters.update_sdp_threshold_max,
				//		Hpp_minus_Hmixed, grad_mixed, find_zeros, external_step, lag);

				El::Matrix<El::BigFloat> search_direction(n_external_parameters, 1);
				for (int i = 0; i < n_external_parameters; i++)
				{
					search_direction(i, 0) = dynamical_parameters.search_direction[i];
				}


				//if (find_zeros && update_sdp)
				if (update_sdp)
				{
					if (compute_ext_step_only_once == true && external_step_specified_Q == true &&
						(recompute_ext_during_re_aiming == false || lowest_mu_Q == false)
						)
					{
						// hess_BFGS will not be updated, but grad_BFGS will be updated
						external_step = external_step_save;
					}
					else
					{
						El::BigFloat lambda;
						El::Matrix<El::BigFloat> external_step_plus, external_step_minus;

						compute_ellipse_boundary(hess_BFGS, grad_corrected, lag_shifted, search_direction,
							external_step_plus, external_step_minus, lambda);
						external_step = external_step_plus;

						if (El::mpi::Rank() == 0)
						{
							if (use_Lpu_mu_correction)
								std::cout << "\nLpu corrected : \n";

							std::cout << "BFGS hess =\n";
							print_matrix(hess_BFGS);

							std::cout << "grad_corrected = ";
							print_vector(grad_corrected);
							std::cout << "\n";

							std::cout << "grad_p = ";
							print_vector(grad_p);
							std::cout << "\n";

							std::cout << "grad_mixed = ";
							print_vector(grad_mixed);
							std::cout << "\n";

							std::cout << " lambda=" << lambda << "\n";
							if (lambda < 0)
								std::cout << "ext-step_zero=";
							else
								std::cout << "ext-step_plus=";

							print_vector(external_step_plus);
							std::cout << "\n" << std::flush;
						}

						external_step_save = external_step;
						external_step_specified_Q = true;
					}

					external_step_size = El::Nrm2(external_step);
				}

			}

			if (beta > 10 && dynamical_parameters.updateSDP_dualityGapThreshold <= 0)
			{
				if (El::mpi::Rank() == 0)std::cout << "beta>10 && (no dualityGapThreshold) : update_sdp=false \n";
				update_sdp = false;
			}

			if (external_step_size < 0 //dynamical_parameters.update_sdp_threshold_min
					//|| (external_step_size > dynamical_parameters.update_sdp_threshold_max && dynamical_parameters.total_iterations == 0)) 
					//|| mu > 1e-10
				//|| (dual_objective) < -0.5 
				)
			{
				if (El::mpi::Rank() == 0)std::cout << "external_step_size < 0 : update_sdp=false \n";
				update_sdp = false;
			}
		}

		//std::cout << "rank=" << El::mpi::Rank() << " pid=" << getpid() << " after if (update_sdp) \n" << std::flush;

		update_dXdY(update_sdp,
			dynamical_parameters, block_info, sdp, grid, X_cholesky, Y_cholesky, timers,
			internal_dx, internal_dy, dx, dy, dX, dY, R,
			delta_lambda, external_step, Delta_xy, primal_step_length, dual_step_length, step_length_reduction);


		if (is_primal_and_dual_feasible)
		{
			primal_step_length = El::Min(primal_step_length, dual_step_length);
			dual_step_length = primal_step_length;
		}

		//El::BigFloat cutoff1 = El::BigFloat(1e-50);
		El::BigFloat cutoff1 = 0;

		if (El::mpi::Rank() == 0)
                {
                        std::cout << "actual beta : " << beta << "\n" << std::flush;
                        if (update_sdp == true)
                        {
                                std::cout << "P/D-step-len-ext = " 
                                          << primal_step_length 
                                          << " , " 
                                          << dual_step_length << "\n" 
                                          << std::flush;
                        }
                }

		if (dynamical_parameters.total_iterations == 0
                   || update_sdp == false
                   || El::Max(primal_step_length, dual_step_length) > PDsteplen_threshold
                   )
                {
                        break;
                }

		if (update_sdp)
		{
			if (stallrecovery_phase == 1 ||
				((El::Max(primal_step_length, dual_step_length) > step_min_threshold &&
				(primal_step_length + dual_step_length) > (primal_step_length_best + dual_step_length_best)
					))
				)
			{
				// save best state
				save_load_beta_scan_best_state(
					dx_best, dy_best, dX_best, dY_best, primal_step_length_best, dual_step_length_best,
					beta_best, hess_BFGS_best, grad_mixed_best,
					dx, dy, dX, dY, primal_step_length, dual_step_length, beta, hess_BFGS, grad_mixed);
			}
		}


		if (El::Max(primal_step_length, dual_step_length) > step_max_threshold) break;

		if ((El::Max(primal_step_length_best, dual_step_length_best) > step_min_threshold || max_climbing <= 0) && beta == El::BigFloat(1))
		{
			// load best state
			save_load_beta_scan_best_state(
				dx, dy, dX, dY, primal_step_length, dual_step_length, beta, hess_BFGS, grad_mixed,
				dx_best, dy_best, dX_best, dY_best, primal_step_length_best, dual_step_length_best,
				beta_best, hess_BFGS_best, grad_mixed_best);

			if (El::mpi::Rank() == 0)std::cout << "reset to best choice with beta = " << beta << "\n" << std::flush;
			break;
		}

		if (El::Max(primal_step_length, dual_step_length) > cutoff1 && beta < 1)
		{
			beta = beta + El::BigFloat(beta_scan_step);
			if (beta > 1)beta = 1;
		}
		// even if beta=1, P/D-step-len-ext still < 0.1 : we re-run internal step with beta=2
		else if (El::Max(primal_step_length, dual_step_length) > cutoff1 && beta == El::BigFloat(1))
		{
			beta = 2;
			max_climbing--;
			lowest_mu_Q = false;
			update_sdp = false;
		}
		// no matter it's in external mode or internal mode, if P/D-step-len is small and beta>1,
		// we accept the internal step and set beta=2 for next step
		else if (El::Max(primal_step_length, dual_step_length) > cutoff1 && beta > El::BigFloat(1))
		{
			intext_mode = 1;
			update_sdp = false;

			update_dXdY(update_sdp,
				dynamical_parameters, block_info, sdp, grid, X_cholesky, Y_cholesky, timers,
				internal_dx, internal_dy, dx, dy, dX, dY, R,
				delta_lambda, external_step, Delta_xy, primal_step_length, dual_step_length, step_length_reduction);

			break;
		}

		if (beta >= El::BigFloat(1e50))break;
		//if (beta >= El::BigFloat(1))break;

		if (El::Max(primal_step_length, dual_step_length) <= cutoff1)
		{
			if (El::mpi::Rank() == 0)std::cout << "max_step.beta=" << El::Max(primal_step_length, dual_step_length)*beta
				<< " max_step=" << El::Max(primal_step_length, dual_step_length)
				<< ". reset beta to 1e50. "
				<< "\n" << std::flush;
			beta = El::BigFloat(1e50);
		}

	}
	external_step *= dual_step_length;
	external_step_size = El::Nrm2(external_step);

	if (centering_mode && El::Min(primal_step_length, dual_step_length) < 0.05)
	{
		intext_mode = 1;
	}

	// update grad_BFGS
	if (update_sdp)
	{
		if (compare_BFGS_gradient_at_same_mu)
		{
			extrapolate_gradient_as_target_mu(primal_step_length, primal_step_length, grad_p, grad_mixed, grad_BFGS);
		}
		else if (use_gradp_for_BFGS_update)
		{
			grad_BFGS = grad_p;
		}
		else
		{
			grad_BFGS = grad_mixed;
			grad_BFGS += grad_p;
		}
	}

	// When we change our decision to update sdp based on the final scaled external_step_size,
	// we redo some computations using internal_dx, internal_dy instead of dx, dy. 
	// For now we reset external_step to be zero but not the external_step_size
	// so that we can do sanity checks in the log file where external_step_size where printed out in each iteration.
	// Else catches 1) previously decided update_sdp = false where dx, dy = internal_dx, internal_dy;
	//              2) update_sdp = true still holds.   
	if (update_sdp
		&& (external_step_size > dynamical_parameters.update_sdp_threshold_max)
		)
		//|| external_step_size < dynamical_parameters.update_sdp_threshold_min))
	{
		if (El::mpi::Rank() == 0)std::cout << "external_step_size > dynamical_parameters.update_sdp_threshold_max : update_sdp=false \n";
		update_sdp = false;
		El::Zero(external_step);
		//external_step_size = 0 ; 
		constraint_matrix_weighted_sum(block_info, sdp, internal_dx, dX);
		dX += primal_residues;

		multiply(dX, Y, dY);
		dY -= R;
		cholesky_solve(X_cholesky, dY);
		dY.symmetrize();
		dY *= El::BigFloat(-1);

		primal_step_length
			= step_length(X_cholesky, dX, step_length_reduction,
				"run.step.stepLength(XCholesky)", timers);
		dual_step_length
			= step_length(Y_cholesky, dY, step_length_reduction,
				"run.step.stepLength(YCholesky)", timers);
		if (is_primal_and_dual_feasible)
		{
			primal_step_length = El::Min(primal_step_length, dual_step_length);
			dual_step_length = primal_step_length;
		}
		for (size_t block = 0; block < x.blocks.size(); ++block)
		{
			El::Axpy(primal_step_length, internal_dx.blocks[block], x.blocks[block]);
		}
		// Update the dual point (y, Y) += dualStepLength*(dy, dY)
		for (size_t block = 0; block < dy.blocks.size(); ++block)
		{
			El::Axpy(dual_step_length, internal_dy.blocks[block], y.blocks[block]);
		}
	}
	else
	{
		for (size_t block = 0; block < x.blocks.size(); ++block)
		{
			El::Axpy(primal_step_length, dx.blocks[block], x.blocks[block]);
		}
		// Update the dual point (y, Y) += dualStepLength*(dy, dY)
		for (size_t block = 0; block < dy.blocks.size(); ++block)
		{
			El::Axpy(dual_step_length, dy.blocks[block], y.blocks[block]);
		}
	}

	predictor_timer.stop();

	dX *= primal_step_length;
	dY *= dual_step_length;
	X += dX;
	Y += dY;
	delta_lambda *= dual_step_length;
	lag_multiplier_lambda += delta_lambda;


	step_timer.stop();
}  
