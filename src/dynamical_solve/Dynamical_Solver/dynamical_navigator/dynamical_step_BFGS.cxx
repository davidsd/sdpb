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

//El::BigFloat compute_lag(const El::BigFloat mu, const Block_Diagonal_Matrix &X_cholesky, const Dynamical_Solver &solver);

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

template <typename T>
void print_stdvector(const std::vector<T> & vec)
{
	if (El::mpi::Rank() == 0)
	{
		for (auto it = vec.begin(); it != vec.end(); it++) std::cout << *it << " ";
	}
}


template <typename T>
void print_stdpair(const std::vector<T> & vec)
{
	if (El::mpi::Rank() == 0)
	{
		for (auto it = vec.begin(); it != vec.end(); it++) std::cout << *it << " ";
	}
}

template <typename T1, typename T2>
std::ostream& operator<<(std::ostream& stream,
	const std::pair<T1,T2>& p) {
	stream << "(" << p.first << "," << p.second << ")";
	return stream;
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

	El::Matrix<El::BigFloat> & grad_withoutlog,
	El::Matrix<El::BigFloat> & grad_withlog,

	std::vector<std::pair<Block_Vector, Block_Vector>> & H_xp,
	std::vector<std::pair<Block_Vector, Block_Vector>> & Delta_xy,
        std::vector<SDP> & eplus_d_sdp);


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

void update_grad_p(const Block_Vector & x,
                   const Block_Vector & y,
                   const std::vector<SDP> & eplus_delta_sdp,
                   const SDP & sdp,
                   const int & n_external_parameters,
                   El::Matrix<El::BigFloat> & grad_p);

void read_prev_grad_step_hess(const Dynamical_Solver_Parameters &dynamical_parameters,
        El::Matrix<El::BigFloat> & prev_grad,
        El::Matrix<El::BigFloat> & prev_step,
        El::Matrix<El::BigFloat> & prev_BFGS,
        El::Matrix<El::BigFloat> & prev_BFGS_pp );


void compute_grad_p_grad_mixed_Hpp_Hmixed(const Dynamical_Solver_Parameters &dynamical_parameters,
	const El::Matrix<El::BigFloat> & eplus, const El::Matrix<El::BigFloat> & eminus, const El::Matrix<El::BigFloat> & esum,
	const std::vector<std::pair<Block_Vector, Block_Vector>> & H_xp, const std::vector<std::pair<Block_Vector, Block_Vector>> & Delta_xy,
	const Block_Vector & internal_dx, const Block_Vector & internal_dy,
	const El::Matrix<El::BigFloat> & grad_withoutlog,
	const El::Matrix<El::BigFloat> & grad_withlog,
	El::Matrix<El::BigFloat> & grad_p, El::Matrix<El::BigFloat> & grad_mixed, El::Matrix<El::BigFloat> & grad_corrected,
	El::Matrix<El::BigFloat> & Lpu, El::BigFloat & mu,
	El::Matrix<El::BigFloat> & hess_pp, El::Matrix<El::BigFloat> & hess_mixed, El::Matrix<El::BigFloat> & hess_Exact);


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


void write_solver_state(const Dynamical_Solver_Parameters &dynamical_parameters,
	const Block_Info &block_info,
	const Block_Diagonal_Matrix &schur_complement_cholesky,
	const Block_Matrix &schur_off_diagonal,
	const El::DistMatrix<El::BigFloat> &Q);

void compute_R_error(const std::size_t &total_psd_rows, const Block_Diagonal_Matrix &X, const Block_Diagonal_Matrix &Y, El::BigFloat & R_error, Timers &timers);


void compute_corrector_R(
	const Block_Info &block_info, const std::size_t &total_psd_rows,

	const Block_Vector &x_const, const Block_Vector &dx_const,
	const Block_Vector &y_const, const Block_Vector &dy_const,
	const Block_Diagonal_Matrix &X_const, const Block_Diagonal_Matrix &dX_const,
	const Block_Diagonal_Matrix &Y_const, const Block_Diagonal_Matrix &dY_const,

	const El::BigFloat &primal_step_length, const El::BigFloat &dual_step_length,

	Block_Diagonal_Matrix &R,
	El::BigFloat &R_error,
	El::BigFloat &mu,

	Timers &timers);

void compute_R_error(const std::size_t &total_psd_rows, const Block_Diagonal_Matrix &X, const Block_Diagonal_Matrix &Y,
	El::BigFloat & R_error, El::BigFloat & mu, Timers &timers);


template <class T1, class T2, class Pred = std::less<T2> >
struct comp_pair_second {
	bool operator()(const std::pair<T1, T2>&left, const std::pair<T1, T2>&right) {
		Pred p;
		return p(left.second, right.second);
	}
};
template <class T1, class T2, class Pred = std::less<T1> >
struct comp_pair_first {
	bool operator()(const std::pair<T1, T2>&left, const std::pair<T1, T2>&right) {
		Pred p;
		return p(left.first, right.first);
	}
};

bool CubicApprox_LineSearch_RecommendPoint(std::vector<std::pair<El::BigFloat, El::BigFloat>> & history,
	const El::BigFloat & beta_scan_end, const El::BigFloat & beta_scan_step, El::BigFloat & new_x);

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


bool update_hess_only_positive = true;

bool use_Lpu_mu_correction = true;

double dgap_near_optimal_threshold = 1e-20;

bool use_gradp_for_BFGS_update = true;

double rescale_initial_hess = 1;

//double BFGS_partial_update_reduction = -1; // if positive, hessian will be updated always positive based on the partial update logic.

bool compare_BFGS_gradient_at_same_mu = false;
bool BFGS_approximate_Hpp_only = true;
//double beta_for_mu_logdetX = 0; // if negative, it will use the beta from beta scan
//double beta_for_mu_logdetX = -1;

/* // SN 2022/12/04 : I moved those parameter to Dynamical_Solver_Parameters

double beta_climbing = 2;
double mu_upper_limit = 0.0001;

double step_min_threshold = 0.1;
double step_max_threshold = 0.6;
double navigator_shift = 0.2;

double beta_scan_begin = -1;
double beta_scan_end = 1.01;
double beta_scan_step = 0.1;
*/

int max_climbing = 1;

El::BigFloat predictor_centering_parameter_V3(const Solver_Parameters &parameters,
	const El::BigFloat & mu, const El::BigFloat & mu_near_optimal_threshold)
{
	return mu < mu_near_optimal_threshold ? parameters.feasible_centering_parameter
		: parameters.infeasible_centering_parameter;
}


El::BigFloat compute_beta_for_finite_dGap(const El::BigFloat & current_dGap, const El::BigFloat & dGap_target)
{
	El::BigFloat beta = dGap_target / current_dGap;
	if (beta < 0.1) beta = El::BigFloat(0.1);
	return beta;
}


bool external_corrector_jump_to_newSDP = true;

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
        if (El::mpi::Rank() == 0) std::cout << "start of step" << "\n";
	auto &step_timer(timers.add_and_start("run.step"));
	step_timer.stop();

	// Internal_step: compute dx and dy for the central sdp as in compute_search_direction()      
	//                - H^-1_xx Del_x L_mu in Eq (12) and Eq(13)
	//                Notice that the negative sign has been accounted. 
	Block_Vector internal_dx(x), internal_dy(y);
	Block_Vector dx(internal_dx), dy(internal_dy);
	Block_Diagonal_Matrix dX(X), dY(Y);
	Block_Vector grad_x(internal_dx), grad_y(internal_dy);
	Block_Diagonal_Matrix R(X);


	// compute Schur complement
	Block_Diagonal_Matrix schur_complement_cholesky(
		block_info.schur_block_sizes(), block_info.block_indices,
		block_info.num_points.size(), grid);
        if (El::mpi::Rank() == 0) std::cout << "FINISH schur_complement_cholesky" << "\n";
	Block_Matrix schur_off_diagonal(sdp.free_var_matrix);
        if (El::mpi::Rank() == 0) std::cout << "schur_off_diagonal" << "\n";
	El::DistMatrix<El::BigFloat> Q(sdp.dual_objective_b.Height(),
		sdp.dual_objective_b.Height());
        if (El::mpi::Rank() == 0) std::cout << "Initialize Q" << "\n";
	if (!(dynamical_parameters.external_corrector_Q && external_corrector_jump_to_newSDP))
	{
                if (El::mpi::Rank() == 0) std::cout << "after external_corrector_Q" << "\n";
		initialize_schur_complement_solver(block_info, sdp, A_X_inv, A_Y, grid,
			schur_complement_cholesky, schur_off_diagonal, Q, timers);
	}
        if (El::mpi::Rank() == 0) std::cout << "after initialize_schur_complement_solver" << "\n";
	auto &frobenius_timer(
		timers.add_and_start("run.step.frobenius_product_symmetric"));
	El::BigFloat dGap_current = frobenius_product_symmetric(X, Y);
	frobenius_timer.stop();

	mu = dGap_current / total_psd_rows;

	if (mu > dynamical_parameters.solver_parameters.max_complementarity)
	{
		terminate_now = true;
		return;
	}

	El::BigFloat beta_scan_begin_El = dynamical_parameters.beta_scan_min;
	El::BigFloat beta_internal = predictor_centering_parameter_V3(dynamical_parameters.solver_parameters, dGap_current, dgap_near_optimal_threshold);
	if (beta_scan_begin_El < 0)beta_scan_begin_El = beta_internal;
	El::BigFloat beta_scan_end_El = dynamical_parameters.beta_scan_max;
	El::BigFloat beta_scan_step_El = dynamical_parameters.beta_scan_step;
	beta = beta_scan_begin_El;

	El::BigFloat step_length_reduction = dynamical_parameters.solver_parameters.step_length_reduction;
	if (dGap_current < dgap_near_optimal_threshold)
	{
		//step_length_reduction = 0.9;
		beta = 0.1;
	}

	/**/
	if (El::mpi::Rank() == 0) std::cout
		<< " duality_gap=" << duality_gap
		<< " dGap_current=" << dGap_current
		<< " updateSDP_dualityGapThreshold=" << dynamical_parameters.updateSDP_dualityGapThreshold
		<< " bool1=" << (dynamical_parameters.updateSDP_dualityGapThreshold > 0 && duality_gap > dynamical_parameters.updateSDP_dualityGapThreshold)
		<< " bool2=" << (duality_gap == El::BigFloat(0))
		<< "\n";


	// ordinary sdpb run until dualityGap < dualityGapThreshold
	// this should only happend in the 1st call
	if ((dynamical_parameters.total_iterations == 0 && dGap_current > 0 &&
		dynamical_parameters.updateSDP_dualityGapThreshold > 0 && dGap_current > dynamical_parameters.updateSDP_dualityGapThreshold) ||
		duality_gap == El::BigFloat(0))
	{
		update_sdp = false;
		return internal_step(dynamical_parameters, total_psd_rows, block_info, sdp, grid,
			X_cholesky, Y_cholesky, timers,
			schur_complement_cholesky, schur_off_diagonal, Q,
			dx, dy, dX, dY, R, grad_x, grad_y,
			primal_residue_p, mu, is_primal_and_dual_feasible,
			beta_internal, primal_step_length, dual_step_length, step_length_reduction);
	}


	if (El::mpi::Rank() == 0)
		std::cout << "R_error : " << R_error << " \n";

	if (dynamical_parameters.external_corrector_Q && external_corrector_jump_to_newSDP)
	{
		update_sdp = false;
		external_corrector_jump_to_newSDP = false;
		return external_corrector_run(dynamical_parameters, sdp,
			block_info, grid,
			dx, dy, dX, dY,
			primal_step_length, dual_step_length, step_length_reduction,
			timers);
	}

	if (El::mpi::Rank() == 0) std::cout
		<< " finite_dGap_target=" << dynamical_parameters.finite_dGap_target
		<< " abs(current-target)=" << El::Abs(dynamical_parameters.finite_dGap_target - dGap_current)
		<< " updateSDP_dualityGapThreshold=" << dynamical_parameters.updateSDP_dualityGapThreshold
		<< " bool1=" << (El::Abs(dynamical_parameters.finite_dGap_target - dGap_current) > dynamical_parameters.centeringRThreshold)
		<< " bool2=" << (dynamical_parameters.centeringRThreshold > 0 && R_error > dynamical_parameters.centeringRThreshold)
		<< "\n";

        //std::vector<SDP> eplus_delta_sdp;
        //read_sdp_grid(dynamical_parameters, block_info, sdp, grid, timers,
        //        schur_complement_cholesky, schur_off_diagonal, Q, x, y, X_cholesky, n_external_parameters,
        //        eplus, eminus, esum, Lpu, grad_withoutlog, grad_withlog, H_xp, Delta_xy,eplus_delta_sdp);

	// move the solver to finite_dGap_target
	if (dynamical_parameters.finite_dGap_target > 0 &&
		(El::Abs(dynamical_parameters.finite_dGap_target - dGap_current) > dynamical_parameters.centeringRThreshold ||
		(dynamical_parameters.centeringRThreshold > 0 && R_error > dynamical_parameters.centeringRThreshold))
		)
	{
		beta = compute_beta_for_finite_dGap(dGap_current, dynamical_parameters.finite_dGap_target);
		update_sdp = false;
		if (El::Min(p_step, d_step) < 0.1)beta = El::Min(1 / (El::Min(p_step, d_step) * 10), El::BigFloat(1e10));

		if (El::mpi::Rank() == 0)std::cout << "run centering steps with beta=" << beta
			<< " duality Gap target=" << dynamical_parameters.finite_dGap_target
			<< " abs(target-current)=" << El::Abs(dynamical_parameters.finite_dGap_target - dGap_current)
			<< " prev_p_step=" << p_step << " prev_d_step=" << d_step
			<< "\n";

		return internal_step(dynamical_parameters, total_psd_rows, block_info, sdp, grid,
			X_cholesky, Y_cholesky, timers,
			schur_complement_cholesky, schur_off_diagonal, Q,
			dx, dy, dX, dY, R, grad_x, grad_y,
			primal_residue_p, mu, is_primal_and_dual_feasible,
			beta, primal_step_length, dual_step_length, step_length_reduction);
	}


	if (dynamical_parameters.finite_dGap_target > 0 && dynamical_parameters.returnCheckpointOnLCP && 
		dynamical_parameters.centeringRThreshold > 0 && R_error < dynamical_parameters.centeringRThreshold &&
		El::Abs(dynamical_parameters.finite_dGap_target - dGap_current) < dynamical_parameters.centeringRThreshold )
	{
		update_sdp = true;

		beta = compute_beta_for_finite_dGap(dGap_current, dynamical_parameters.finite_dGap_target);
		update_sdp = false;
		if (El::Min(p_step, d_step) < 0.1)beta = El::Min(1 / (El::Min(p_step, d_step) * 10), El::BigFloat(1e10));

		return internal_step(dynamical_parameters, total_psd_rows, block_info, sdp, grid,
			X_cholesky, Y_cholesky, timers,
			schur_complement_cholesky, schur_off_diagonal, Q,
			dx, dy, dX, dY, R, grad_x, grad_y,
			primal_residue_p, mu, is_primal_and_dual_feasible,
			beta, primal_step_length, dual_step_length, step_length_reduction);
	}


	// if R_error is not small, we do centering steps
	if (dynamical_parameters.updateSDP_dualityGapThreshold <= 0 &&
		dynamical_parameters.centeringRThreshold > 0 && R_error > dynamical_parameters.centeringRThreshold)
	{
		if (El::mpi::Rank() == 0)std::cout << "run centering steps with beta=1 \n";
		beta = 1;
		update_sdp = false;

		return internal_step(dynamical_parameters, total_psd_rows, block_info, sdp, grid,
			X_cholesky, Y_cholesky, timers,
			schur_complement_cholesky, schur_off_diagonal, Q,
			dx, dy, dX, dY, R, grad_x, grad_y,
			primal_residue_p, mu, is_primal_and_dual_feasible,
			beta, primal_step_length, dual_step_length, step_length_reduction);
	}

	// from here, the solver is on the central path. We must try to hop.

	Block_Diagonal_Matrix dX_best(X), dY_best(Y);
	Block_Vector dx_best(x), dy_best(y);
	El::BigFloat primal_step_length_best, dual_step_length_best, beta_best;
	El::BigFloat primal_step_length_prev, dual_step_length_prev, beta_prev;
	El::Matrix<El::BigFloat> hess_BFGS_best, grad_mixed_best;

	int n_external_parameters = dynamical_parameters.n_external_parameters;
	std::vector<std::pair<Block_Vector, Block_Vector>> H_xp, Delta_xy;
	El::Matrix<El::BigFloat> eplus(n_external_parameters, 1), eminus(n_external_parameters, 1), esum(n_external_parameters, n_external_parameters);
	El::Matrix<El::BigFloat> Lpu(n_external_parameters, 1);

	El::Matrix<El::BigFloat> grad_p(n_external_parameters, 1); // Del_p L
	El::Matrix<El::BigFloat> hess_pp(n_external_parameters, n_external_parameters); //H_pp

	El::Matrix<El::BigFloat> grad_corrected(n_external_parameters, 1);

	El::Matrix<El::BigFloat> hess_mixed(n_external_parameters, n_external_parameters); //H_px H^-1_xx H_xp in Eq(13).
	//El::Matrix<El::BigFloat> grad_mixed(n_external_parameters, 1); //H_px (-internal_dx_dy) = H_px (H^-1_xx Del_x L_mu) in Eq (13).  

	El::Matrix<El::BigFloat> prev_grad(n_external_parameters, 1), prev_step(n_external_parameters, 1), prev_BFGS(n_external_parameters, n_external_parameters) , prev_BFGS_pp(n_external_parameters, n_external_parameters);
	El::Zeros(hess_Exact, n_external_parameters, n_external_parameters);
  El::Zeros(grad_mixed, n_external_parameters, 1);
  

	read_prev_grad_step_hess(dynamical_parameters, prev_grad, prev_step, prev_BFGS, prev_BFGS_pp);



  std::vector<SDP> eplus_delta_sdp;
	read_sdp_grid(dynamical_parameters, block_info, sdp, grid, timers,
		schur_complement_cholesky, schur_off_diagonal, Q, x, y, X_cholesky, n_external_parameters,
		eplus, eminus, esum, Lpu, grad_withoutlog, grad_withlog, H_xp, Delta_xy,eplus_delta_sdp);


	if (El::mpi::Rank() == 0)std::cout << "eplus[0]=" << eplus(0) 
		<< " grad_withoutlog=" << grad_withoutlog(0)
		<< " Lpu(0)=" << Lpu(0) 
		<< "\n";

	// beta scan
	for (beta = beta_scan_begin_El; beta <= beta_scan_end_El; beta += beta_scan_step_El)
	{
		internal_predictor_direction_dxdydXdY(block_info, sdp, *this, schur_complement_cholesky,
			schur_off_diagonal, X_cholesky, beta,
			mu, primal_residue_p, Q, grad_x, grad_y, internal_dx, internal_dy, dX, dY, R);
		internal_corrector_direction(block_info, sdp, *this, schur_complement_cholesky,
			schur_off_diagonal, X_cholesky, beta,
			mu, primal_residue_p, Q, grad_x, grad_y, internal_dx, internal_dy, dX, dY, R);
                update_Lpu(block_info, sdp, internal_dx, X_cholesky);
		// compute various variables according to the formula
                // grad_corrected = grad_p (grad_withoutlog or grad_withlog) 
                //                  - mu Lpu  (option 1)
                //                  + grad_mixed (option 2) 
                // Lpu = d(Log det X)/dp, To be chekced: grad_mixed - (- mu Lpu) ~ O(R) 
                // options decided by use_Lpu_mu_correction (false by default) 
                // grad_withoutlog = dc.x + c.dx 
                // grad_withlog    = db.y + dc.x + x.dB.y
                // decided by dynamical_parameters.gradientWithLogDetX
		compute_grad_p_grad_mixed_Hpp_Hmixed(dynamical_parameters,
			eplus, eminus, esum, H_xp, Delta_xy, internal_dx, internal_dy,
			grad_withoutlog, grad_withlog, grad_p, grad_mixed, grad_corrected,
			Lpu, mu, hess_pp, hess_mixed, hess_Exact);
                //El::Zero(hess_BFGS_pp); 
                hess_BFGS_pp = hess_Exact; 
                if (dynamical_parameters.use_exact_hessian){
                   hess_BFGS_pp = hess_pp;  
                 }    
		// decide hess_BFGS
                // grad_diff = grad_p (option 1)
                //           = grad_corrected (option 2)
                // decided by use_gradp_for_BFGS_update (true by default)
                // hess_BFGS (grad_p) ~ hess_BFGS (grad_corrected)   
                // H_pp - H_px H_xx^-1 H_xp ~ 
                // grad_diff = grad_p(i+1) + 0 - (grad_p(i) + step_length grad_mixed(i)) 
                // Assuming R = 0, grad_diff = grad_p(i+1) + 0 - (grad_p(i) + step_length * mu Lpu ) 
		strategy_hess_BFGS(dynamical_parameters, n_external_parameters,
			lowest_mu_Q, grad_p, grad_mixed, grad_corrected,
			Lpu, mu, hess_pp, hess_mixed, hess_Exact,
			prev_BFGS,prev_BFGS_pp, prev_step, prev_grad, hess_BFGS_lowest_mu);

		// compute external step
		if (dynamical_parameters.find_boundary && dynamical_parameters.total_iterations > 0)
		{
			strategy_findboundary_extstep(block_info, X_cholesky, total_psd_rows, mu, beta,
				n_external_parameters, dynamical_parameters, lowest_mu_Q,
				grad_corrected, grad_p, grad_mixed,
				external_step, external_step_save, external_step_specified_Q);
		}
		else
		{
			findMinimumQ = true;
			// this is not used in minimization mode. we should compute lag_shifted only in "debug" mode
			lag_shifted = finite_mu_navigator(block_info, X_cholesky, total_psd_rows, mu, beta, dynamical_parameters);

			external_step = grad_corrected;
			external_step *= (-1);
			El::LinearSolve(hess_BFGS, external_step);

			if (El::mpi::Rank() == 0)
			{
				std::cout << "Hessian Mixed =\n";
				print_matrix(hess_mixed);
				std::cout << "ext-step=";
				print_vector(external_step);
                                std::cout << "grad_corrected= ";
                                print_vector(grad_corrected);
			}
		}

		// compute dx,dy,dX,dY for external_step
		compute_external_dxdydXdY(is_primal_and_dual_feasible,
			dynamical_parameters, block_info, sdp, grid, X_cholesky, Y_cholesky, timers,
			internal_dx, internal_dy, dx, dy, dX, dY, R,
			external_step, Delta_xy, primal_step_length, dual_step_length, step_length_reduction);

		// save the best beta
		if (beta == beta_scan_begin_El ||
			El::Max(primal_step_length, dual_step_length) > El::Max(primal_step_length_best, dual_step_length_best)
			)
		{
			// save best state
			save_load_beta_scan_best_state(
				dx_best, dy_best, dX_best, dY_best, primal_step_length_best, dual_step_length_best,
				beta_best, hess_BFGS_best, grad_mixed_best,
				dx, dy, dX, dY, primal_step_length, dual_step_length, beta, hess_BFGS, grad_mixed);
		}

		if (El::mpi::Rank() == 0)
		{
			std::cout << "scan beta : " << beta << "\n" << std::flush;
			std::cout << "P/D-step-len-ext = " << primal_step_length << " , " << dual_step_length << "\n" << std::flush;
		}

		// there is no special reason we don't want re-aiming for 1st call.
		// But I just prefer to seperate them (usually 2nd call has some climing steps)
		if (dynamical_parameters.total_iterations == 0)
		{
			update_sdp = true;
			execute_step(dx, dy, dX, dY, primal_step_length, dual_step_length);
			external_step *= dual_step_length;
			external_step_size = El::Nrm2(external_step);
			strategy_update_grad_BFGS(primal_step_length, dual_step_length, grad_p, grad_mixed, grad_BFGS);
			return;
		}

		// if step_length > step_max_threshold, break
		if (El::Max(primal_step_length, dual_step_length) > dynamical_parameters.step_max_threshold) break;
	}

	// load best state
	save_load_beta_scan_best_state(
		dx, dy, dX, dY, primal_step_length, dual_step_length, beta, hess_BFGS, grad_mixed,
		dx_best, dy_best, dX_best, dY_best, primal_step_length_best, dual_step_length_best,
		beta_best, hess_BFGS_best, grad_mixed_best);

	if (El::mpi::Rank() == 0)std::cout << "best scan result : best beta = " << beta << " step-len : ("
		<< primal_step_length << ", " << dual_step_length << ")\n" << std::flush;

	if (El::mpi::Rank() == 0)std::cout << "dualityGap_upper_limit=" << dynamical_parameters.dualityGap_upper_limit
		<< " 2 * mu*total_psd_rows=" << 2 * mu*total_psd_rows
		<< "\n" << std::flush;

	// check beta scan result
	if (El::Max(primal_step_length, dual_step_length) > dynamical_parameters.step_min_threshold || max_climbing <= 0 ||
		(dynamical_parameters.dualityGap_upper_limit > 0 && 2 * mu*total_psd_rows > dynamical_parameters.dualityGap_upper_limit)
		)
	{
		// hopping step

		if (dynamical_parameters.external_corrector_Q)
		{
			save_solver_state(dynamical_parameters, block_info,
				dx, dy, dX, dY,
				schur_complement_cholesky, schur_off_diagonal, Q);
		}

		El::BigFloat R_err, R_err2, coit_mu, coit_mu2;
		Block_Diagonal_Matrix R_0(X);
		compute_corrector_R(block_info, total_psd_rows,
			x, dx, y, dy, X, dX, Y, dY,
			1, 1, R_0, R_err, coit_mu, timers);
		compute_R_error(total_psd_rows, X, Y, R_err2, coit_mu2, timers);
		if (El::mpi::Rank() == 0)std::cout << "before hopping step, R_err(X)=" << R_err2
			<< " mu=" << coit_mu2
			<< " R_err(X+dX)=" << R_err
			<< " mu=" << coit_mu << "\n" << std::flush;

		update_sdp = true;
		execute_step(dx, dy, dX, dY, primal_step_length, dual_step_length);
                //update_grad_p(x, y, eplus_delta_sdp, sdp, n_external_parameters, grad_p);
		compute_R_error(total_psd_rows, X, Y, R_err, coit_mu, timers);
		if (El::mpi::Rank() == 0)std::cout << "after hopping step, R_err = " << R_err
			<< " mu=" << coit_mu << "\n" << std::flush;

		external_step *= dual_step_length;
		external_step_size = El::Nrm2(external_step);
		strategy_update_grad_BFGS(primal_step_length, dual_step_length, grad_p, grad_mixed, grad_BFGS);
		return;
	}

	// climbing
	// internal step
	update_sdp = false;

	beta = dynamical_parameters.beta_climbing;
	max_climbing--;
	lowest_mu_Q = false;

	return internal_step(dynamical_parameters, total_psd_rows, block_info, sdp, grid,
		X_cholesky, Y_cholesky, timers,
		schur_complement_cholesky, schur_off_diagonal, Q,
		dx, dy, dX, dY, R, grad_x, grad_y,
		primal_residue_p, mu, is_primal_and_dual_feasible,
		beta, primal_step_length, dual_step_length, step_length_reduction);

}
