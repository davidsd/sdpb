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


void find_zero_step (const El::BigFloat &thresh, const int &max_it, const El::BigFloat &step_size_max,
		El::Matrix<El::BigFloat> &hess, const El::Matrix<El::BigFloat> &grad, bool &find_zeros, 
		El::Matrix<El::BigFloat> &external_step, const El::BigFloat &quadratic_const);





void BFGS_update_hessian(const int n_parameters,
		const El::Matrix<El::BigFloat> &grad_p_diff,
		const El::Matrix<El::BigFloat> &last_it_step,
		El::Matrix<El::BigFloat> &hess_bfgs
		);

void Dynamical_Solver::update_dXdY(bool external_step_Q,

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

	El::BigFloat &delta_lambda,
	El::Matrix<El::BigFloat> & external_step,
	std::vector<std::pair<Block_Vector, Block_Vector>> & Delta_xy,
	El::BigFloat &primal_step_length, El::BigFloat &dual_step_length)
{
	dx = internal_dx;
	dy = internal_dy;
	// Update dx and dy if the external parameter step is small enough 
	// RHS of Eq(12)
	//    - H^-1_xx Del_x L_mu - H^-1_xx H_xp dp 
	// =  internal_dx, internal_dy + Delta_xy . external_step
	if (external_step_Q)
	{
		for (int i = 0; i<dynamical_parameters.n_external_parameters; i++)
		{
			for (size_t block = 0; block < x.blocks.size(); ++block)
			{
				if (dynamical_parameters.find_boundary)
				{
					dx.blocks[block] *= (1.0 + delta_lambda / lag_multiplier_lambda);
				}
				El::Axpy(external_step(i), Delta_xy.at(i).first.blocks[block], dx.blocks[block]);
			}
			for (size_t block = 0; block < dy.blocks.size(); ++block)
			{
				if (dynamical_parameters.find_boundary)
				{
					dy.blocks[block] *= (1.0 + delta_lambda / lag_multiplier_lambda);
				}
				El::Axpy(external_step(i), Delta_xy.at(i).second.blocks[block], dy.blocks[block]);
			}
		}

		if (dynamical_parameters.find_boundary)
		{
			primal_residues *= (1.0 + delta_lambda / lag_multiplier_lambda);
			R *= (1.0 + delta_lambda / lag_multiplier_lambda);
		}
	}
	//Block_Diagonal_Matrix dX(X), dY(Y);
	// dX = PrimalResidues + \sum_p A_p dx[p]
	constraint_matrix_weighted_sum(block_info, sdp, dx, dX);
	dX += primal_residues;

	// dY = Symmetrize(X^{-1} (R - dX Y))
	multiply(dX, Y, dY);
	dY -= R;
	cholesky_solve(X_cholesky, dY);
	dY.symmetrize();
	dY *= El::BigFloat(-1);

	//std::cout << "rank=" << El::mpi::Rank() << " before primal_step_length \n" << std::flush;

	// Compute step-lengths that preserve positive definiteness of X, Y
	primal_step_length
		= step_length(X_cholesky, dX, dynamical_parameters.solver_parameters.step_length_reduction,
			"run.step.stepLength(XCholesky)", timers);

	//std::cout << "rank=" << El::mpi::Rank() << " after primal_step_length \n" << std::flush;

	dual_step_length
		= step_length(Y_cholesky, dY, dynamical_parameters.solver_parameters.step_length_reduction,
			"run.step.stepLength(YCholesky)", timers);
}



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

	int n_external_parameters,

	El::Matrix<El::BigFloat> & eplus,
	El::Matrix<El::BigFloat> & eminus,
	El::Matrix<El::BigFloat> & esum,

	std::vector<std::pair<Block_Vector, Block_Vector>> & H_xp,
	std::vector<std::pair<Block_Vector, Block_Vector>> & Delta_xy)
{
	El::Zero(eplus);
	El::Zero(eminus);
	El::Zero(esum);
	if (dynamical_parameters.new_sdp_path.extension() == ".nsv")
	{
		for (auto &filename : read_file_list(dynamical_parameters.new_sdp_path))
		{
			//Assume that the filename takes the form "plus_i","minus_i" and "sum_i_j", f
			//standing for the change in positive e_i, negative e_i and (e_i + e_j) directions respectively 
			std::string file_name = filename.stem().string();
			std::vector<std::string> directions;
			boost::algorithm::split(directions, file_name, boost::is_any_of("_"));
			SDP new_sdp(filename, block_info, grid), d_sdp(new_sdp);
			Axpy(El::BigFloat(-1), sdp, d_sdp);

			//std::cout << "rank=" << El::mpi::Rank() << " pid=" << getpid() << " before approx_obj \n" << std::flush;

			Approx_Objective approx_obj(block_info, sdp, d_sdp, x, y,
				schur_complement_cholesky,
				schur_off_diagonal, Q);

			//std::cout << "rank=" << El::mpi::Rank() << " after approx_obj \n" << std::flush;

			if (directions[0] == "plus")
			{
				mixed_hess(block_info, d_sdp, x, y, schur_complement_cholesky, schur_off_diagonal, Q, dynamical_parameters.alpha, H_xp, Delta_xy);
				eplus(std::stoi(directions[1])) = approx_obj.d_objective;//+ approx_obj.dd_objective; 
																		 //compute_delta_lag(d_sdp, *this);
			}
			else if (directions[0] == "minus")
			{
				eminus(std::stoi(directions[1])) = approx_obj.d_objective;//+ approx_obj.dd_objective; 
																		  //compute_delta_lag(d_sdp,*this); 
			}
			else if (directions[0] == "sum")
			{
				El::BigFloat tempt = approx_obj.d_objective;//+ approx_obj.dd_objective;
				esum(std::stoi(directions[1]), std::stoi(directions[2])) = tempt;
				esum(std::stoi(directions[2]), std::stoi(directions[1])) = tempt;
				//compute_delta_lag(d_sdp,*this); 
			}
		}
	}
	else
	{
		throw std::invalid_argument("A list of perturbed sdp files are required");
	}
	return;
}

