#include <unistd.h>

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
	for(int i = 0;i<n_parameters;i++){
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


void find_zero_step (const El::BigFloat &thresh, const int &max_it, const El::BigFloat &step_size_max,
		El::Matrix<El::BigFloat> &hess, const El::Matrix<El::BigFloat> &grad, bool &find_zeros, 
		El::Matrix<El::BigFloat> &external_step, const El::BigFloat &quadratic_const);


// Centering parameter \beta_p for the predictor step
El::BigFloat predictor_centering_parameter_V2(const Solver_Parameters &parameters,
		const bool is_primal_dual_feasible)
{
	return is_primal_dual_feasible ? parameters.feasible_centering_parameter
		: parameters.infeasible_centering_parameter;
}



void print_matrix(El::Matrix<El::BigFloat> & matrix)
{
	int dim_width = matrix.Width();
	int dim_height = matrix.Height();

	for (int i = 0; i < dim_height; i++)
	{
		for (int j = 0; j < dim_width; j++)if (El::mpi::Rank() == 0) std::cout << std::setprecision(100) << matrix(i, j) << " ";
		if (El::mpi::Rank() == 0) std::cout << "\n";
	}
	return;
}



void external_grad(const El::Matrix<El::BigFloat> &ePlus,
	const El::Matrix<El::BigFloat> &eMinus,
	const El::BigFloat &alpha,
	El::Matrix<El::BigFloat> &grad)
{
	/*
	grad = ePlus;
	grad -= eMinus;
	grad *= El::BigFloat(1) / (El::BigFloat(2)*alpha);
	*/

	grad = ePlus;
	grad *= El::BigFloat(1) / alpha;
}



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

	int n_external_parameters,

	El::Matrix<El::BigFloat> & eplus,
	El::Matrix<El::BigFloat> & eminus,
	El::Matrix<El::BigFloat> & esum,

	std::vector<std::pair<Block_Vector, Block_Vector>> & H_xp,
	std::vector<std::pair<Block_Vector, Block_Vector>> & Delta_xy);


// A subroutine called by run_dynamical
// The external parameter step is passed by the argument 'external_step'
// Compute external_step using Eq (13)
// Compute (dx, dy, dX, dY) using Eq(12)
// Scale both external_step and (dx, dy, dX, dY) by the step length 

int centering_steps = 5;

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

	bool jump(false);

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
	if(mu > dynamical_parameters.solver_parameters.max_complementarity)
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
	beta  = predictor_centering_parameter_V2(dynamical_parameters.solver_parameters, is_primal_and_dual_feasible);

	//std::cout << "rank=" << El::mpi::Rank() << " pid=" << getpid() << " check R_error=" << R_error << "\n" << std::flush;

	if (dynamical_parameters.updateSDP_dualityGapThreshold <= 0 && 
		dynamical_parameters.centeringRThreshold>0 && R_error > dynamical_parameters.centeringRThreshold)
	{
		if (El::mpi::Rank() == 0)std::cout << "run centering steps with beta=1 \n";
		beta = 1;
		update_sdp = false;
	}

	//if (mu_direction_mode == 1)beta = 1;
	//if (mu_direction_mode == 2)beta = 2;

	//std::cout << "status : intext_mode=" << intext_mode << " beta=" << beta << " update_sdp=" << update_sdp <<  " \n";

	if (intext_mode == 1)
	{
		beta = 2;
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

	//if (El::mpi::Rank() == 0)std::cout << "test 1 \n" << std::flush;

	//std::cout << "rank=" << El::mpi::Rank() << " pid=" << getpid() << " before for, update_sdp=" << update_sdp << "\n" << std::flush;

	int n_external_parameters = dynamical_parameters.n_external_parameters;
	std::vector<std::pair<Block_Vector, Block_Vector>> H_xp, Delta_xy;
	El::Matrix<El::BigFloat> eplus(n_external_parameters, 1), eminus(n_external_parameters, 1), esum(n_external_parameters, n_external_parameters);

	if (update_sdp)
		read_sdp_grid(dynamical_parameters, block_info, sdp, grid, timers,
			schur_complement_cholesky, schur_off_diagonal, Q, x, y, n_external_parameters,
			eplus, eminus, esum,
			H_xp, Delta_xy);

	for (int stallrecovery_phase = 1; stallrecovery_phase <= 20; stallrecovery_phase++)
	{
		//if (El::mpi::Rank() == 0)std::cout << "test 2 \n" << std::flush;

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
			El::Matrix<El::BigFloat> grad_p(n_external_parameters,1); // Del_p L
			El::Matrix<El::BigFloat> hess_pp(n_external_parameters, n_external_parameters); //H_pp

			if (dynamical_parameters.use_exact_hessian)
			{
				external_grad_hessian(eplus, eminus, esum, dynamical_parameters.alpha, grad_p, hess_pp);
			}
			else
			{
				//std::cout << "compute grad: " << '\n' << std::flush;
				external_grad(eplus, eminus, dynamical_parameters.alpha, grad_p);
			}

			El::Matrix<El::BigFloat> hess_mixed(n_external_parameters,n_external_parameters); //H_px H^-1_xx H_xp in Eq(13).
			El::Matrix<El::BigFloat> grad_mixed(n_external_parameters,1); //H_px (-internal_dx_dy) = H_px (H^-1_xx Del_x L_mu) in Eq (13).  

			for (int i = 0; i<n_external_parameters; i++)
			{
				for (int j = 0; j<n_external_parameters; j++)
				{
					// The minus sign compensate the minus sign when calculating Delta_xy in Eq(15)
					hess_mixed(i, j) = -(dot(H_xp.at(i).first, Delta_xy.at(j).first) + dot(H_xp.at(i).second, Delta_xy.at(j).second));
				}
			}

			for (int i=0; i<n_external_parameters; i++)
			{ 
				// The minus sign compensates the minus sign when calculating the internal step 
				grad_mixed(i) = (dot(H_xp.at(i).first,internal_dx) + dot(H_xp.at(i).second,internal_dy));
			}  

			// Eq(13):
			// (H_pp - H_px H^-1_xx H_xp)      dp     =  - Del_p L   + H_px (H^-1_xx Del_x L)
			// (H_pp - (H_px H^-1_xx H_xp))    dp     =  - Del_p L   + H_px . (- internal_dx, - internal_dy)  
			// (hess_pp - hess_mixed)          dp     =  - grad_p    + grad_mixed 

			// H_pp - H_px H^-1_xx H_xp, LHS of Eq(13)
			El::Matrix<El::BigFloat> Hpp_minus_Hmixed(hess_pp);

			if (dynamical_parameters.use_exact_hessian)
			{
				Hpp_minus_Hmixed -= hess_mixed;
				hess_Exact = Hpp_minus_Hmixed;
			}

			// - H_px (H^-1_xx Del_x L_mu) + Del_p L_mu, RHS of Eq(13) 
			grad_mixed += grad_p;
			El::Matrix<El::BigFloat> &grad_corrected = grad_mixed;

			grad_BFGS = grad_mixed;
			El::Zeros(hess_BFGS, n_external_parameters,n_external_parameters);
			hess_BFGS -= hess_mixed;

			//std::cout << "hess_mixed: " << '\n' << std::flush;
			//print_matrix(hess_mixed);
			//std::cout << '\n' << std::flush;

			if ( dynamical_parameters.use_exact_hessian
					//|| (mu < 1e-8 && mu > 1e-9)
			   )
			{
				//Initialize the BFGS hessian for the next
				hess_BFGS = Hpp_minus_Hmixed;

				//Flip the sign of Hessian if determinant is negative 
				typedef El::Base<El::BigFloat> Real;
				El::Matrix<Real> w(n_external_parameters,1);
				El::Zero(w);
				El::Matrix<El::BigFloat> Q(Hpp_minus_Hmixed);
				El::Matrix<El::BigFloat> tempt(Hpp_minus_Hmixed);
				El::Matrix<El::BigFloat> diag(Hpp_minus_Hmixed);
				El::Zero(diag);
				El::HermitianEig(El::LOWER, Hpp_minus_Hmixed, w, Q);//,ctrl);
				for (int i = 0; i < n_external_parameters; i++){
					//std::cout<< "diag: " << w(i) << ", " << std::flush;
					if (w(i) < 0){if (El::mpi::Rank() == 0)std::cout<< "flip the sign of hessian" <<'\n' << std::flush;}
					diag(i,i) = El::Abs(w(i));
				}
				El::Gemm(El::NORMAL, El::TRANSPOSE, El::BigFloat(1), diag, Q, tempt);
				El::Gemm(El::NORMAL, El::NORMAL, El::BigFloat(1),  Q, tempt, Hpp_minus_Hmixed);

				external_step = grad_mixed;
				external_step *= (-1);// we can get rid of grad_mixed here but it might be confusing. 

				El::LinearSolve(hess_Exact, external_step);
			}

			else{
				if (El::mpi::Rank() == 0)std::cout<<"using BFGS: "<<'\n'<<std::flush;

				El::Matrix<El::BigFloat> prev_grad(n_external_parameters,1), prev_step(n_external_parameters,1);
				for (int i=0; i<n_external_parameters; i++)
				{ prev_grad(i,0) = dynamical_parameters.prev_grad[i];
					prev_step(i,0) = dynamical_parameters.prev_step[i];
					if (dynamical_parameters.total_iterations > 0){
						for (int j=0; j<n_external_parameters; j++)
						{
							hess_BFGS(i,j) = dynamical_parameters.hess_BFGS[i * n_external_parameters + j];
						}
					}
				}

				//std::cout << "init BFGS: " << '\n' << std::flush;
				//print_matrix(hess_BFGS);

				if (dynamical_parameters.total_iterations > 0)
				{
					El::Matrix<El::BigFloat> grad_diff; 
					grad_diff = grad_mixed;
					grad_diff -= prev_grad; 
					BFGS_update_hessian(n_external_parameters, grad_diff,prev_step,hess_BFGS);
					//std::cout<<"precision_BFGS_H: "<< hess_BFGS(0,0).Precision() << std::flush; 
				}
				//Flip the sign of Hessian if determinant is negative 
				typedef El::Base<El::BigFloat> Real;
				El::Matrix<Real> w(n_external_parameters,1);
				El::Zero(w);
				El::Matrix<El::BigFloat> Q(hess_BFGS);
				El::Matrix<El::BigFloat> tempt(hess_BFGS);
				El::Matrix<El::BigFloat> diag(hess_BFGS);
				El::Zero(diag);
				El::HermitianEig(El::LOWER, hess_BFGS, w, Q);//,ctrl);
				for (int i = 0; i < n_external_parameters; i++){
					if (w(i) < 0){if (El::mpi::Rank() == 0)std::cout<< "flip the sign of hessian" <<'\n' << std::flush;}
					diag(i,i) = El::Abs(w(i));
				}
				El::Gemm(El::NORMAL, El::TRANSPOSE, El::BigFloat(1), diag, Q, tempt);
				El::Gemm(El::NORMAL, El::NORMAL, El::BigFloat(1),  Q, tempt, hess_BFGS); 

				//std::cout << "updated BFGS: " << '\n' << std::flush;
				//print_matrix(hess_BFGS);

				external_step = grad_mixed; 
				external_step *= (-1);// we can get rid of grad_mixed here but it might be confusing. 
				El::LinearSolve(hess_BFGS, external_step);
				Hpp_minus_Hmixed = hess_BFGS;
			}

			El::Matrix<El::BigFloat> &d_min_p = external_step;
			external_step_size = El::Nrm2(external_step);

			if (dynamical_parameters.find_boundary && El::Abs(duality_gap) < 0.9 && dynamical_parameters.total_iterations > 0)
			{ 
				El::BigFloat lag = compute_lag(beta*mu, X_cholesky, *this);
				El::BigFloat lag0 = compute_lag(0, X_cholesky, *this);
				lag = lag0;

				//if (El::mpi::Rank() == 0)std::cout << "lag  = " << lag << '\n' << std::flush;
				if (El::mpi::Rank() == 0)std::cout << "lag0 = " << lag0 << '\n' << std::flush;

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
					{
						El::BigFloat a, b, c;

						El::Matrix<El::BigFloat> H_pp_d_min_p;
						El::Gemv(El::NORMAL, El::BigFloat(1), hess_pp, d_min_p, H_pp_d_min_p);
						c = El::Dot(d_min_p, grad_mixed) + (dot(grad_x, internal_dx) + dot(grad_y, internal_dy) + lag) + El::Dot(d_min_p, H_pp_d_min_p) / 2;

						El::Matrix<El::BigFloat> Hpp_minus_Hmixed_inv_vp = search_direction;
						El::LinearSolve(Hpp_minus_Hmixed, Hpp_minus_Hmixed_inv_vp);

						El::Matrix<El::BigFloat> Hpp_Hpp_minus_Hmixed_inv_vp;
						El::Gemv(El::NORMAL, El::BigFloat(1), hess_pp, Hpp_minus_Hmixed_inv_vp, Hpp_Hpp_minus_Hmixed_inv_vp);

						a = El::Dot(Hpp_minus_Hmixed_inv_vp, Hpp_Hpp_minus_Hmixed_inv_vp) / 2;

						b = -El::Dot(d_min_p, search_direction) + El::Dot(d_min_p, Hpp_Hpp_minus_Hmixed_inv_vp);

						El::BigFloat b2_minus_4ac = b * b - 4 * a*c;
						El::BigFloat sq = 0;
						if (b2_minus_4ac > 0)sq = El::Sqrt(b2_minus_4ac);
						El::BigFloat lambda_plus = (-b + sq) / (2 * a);
						El::BigFloat lambda_minus = (-b - sq) / (2 * a);

						//std::cout << "a=" << a << " b=" << b << " c=" << c << " b^2-4ac=" << b2_minus_4ac << "\n" << std::flush;

<<<<<<< HEAD
						//std::cout << "lambda_plus=" << lambda_plus << " lambda_minus=" << lambda_minus << "\n" << std::flush;

						El::Matrix<El::BigFloat> external_step_plus = Hpp_minus_Hmixed_inv_vp;
						external_step_plus *= lambda_plus;
						external_step_plus += d_min_p;
						std::cout << "external_step_plus="; for (int i = 0; i < external_step_plus.Height(); i++) std::cout << external_step_plus(i, 0) << " ";
						std::cout << "\n" << std::flush;

						El::Matrix<El::BigFloat> external_step_minus = Hpp_minus_Hmixed_inv_vp;
						external_step_minus *= lambda_minus;
						external_step_minus += d_min_p;
						std::cout << "external_step_minus="; for (int i = 0; i < external_step_minus.Height(); i++) std::cout << external_step_minus(i, 0) << " ";
						std::cout << "\n" << std::flush;
						if (El::mpi::Rank() == 0)std::cout << "new external_step_plus=";
						for (int i = 0; i < external_step_plus.Height(); i++) if (El::mpi::Rank() == 0)std::cout << external_step_plus(i, 0) << " ";
						if (El::mpi::Rank() == 0)std::cout << "\n" << std::flush;

						// in this convention, I believe I should use lambda_plus
						external_step = external_step_plus;

						//std::cout << "external_step="; for (int i = 0; i < external_step.Height(); i++) std::cout << external_step(i, 0) << " ";
						//std::cout << "\n" << std::flush;

						external_step_size = El::Nrm2(external_step);
					}
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

		//if (El::mpi::Rank() == 0)std::cout << "test 5 \n" << std::flush;

		update_dXdY(update_sdp,
			dynamical_parameters, block_info, sdp, grid, X_cholesky, Y_cholesky, timers,
			internal_dx, internal_dy, dx, dy, dX, dY, R,
			delta_lambda, external_step, Delta_xy, primal_step_length, dual_step_length);
			

		//if (El::mpi::Rank() == 0)std::cout << "test 6 \n" << std::flush;

		if(is_primal_and_dual_feasible)
		{
			primal_step_length = El::Min(primal_step_length, dual_step_length);
			dual_step_length = primal_step_length;
		}

		El::BigFloat cutoff1 = El::BigFloat(1e-10);

		if (El::mpi::Rank() == 0)std::cout << "actual beta : " << beta << "\n" << std::flush;

		if (El::mpi::Rank() == 0 && update_sdp == true)
			std::cout << "P/D-step-len-ext = " << primal_step_length << " , " << dual_step_length << "\n" << std::flush;

		if (dynamical_parameters.total_iterations == 0) break;

		if (El::Max(primal_step_length, dual_step_length) > 0.1) break;

		if (El::Max(primal_step_length, dual_step_length) > cutoff1 && beta < 1)
		{
			beta = beta + El::BigFloat(0.1);
			if (beta > 1)beta = 1;
		}
		// even if beta=1, P/D-step-len-ext still < 0.1 : we re-run internal step with beta=2
		else if (El::Max(primal_step_length, dual_step_length) > cutoff1 && beta == El::BigFloat(1))
		{
			beta = 2;
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
				delta_lambda, external_step, Delta_xy, primal_step_length, dual_step_length);

			break;
		}

		if (beta >= El::BigFloat(1e50))break;
		//if (beta >= El::BigFloat(1))break;
  
		if(El::Max(primal_step_length, dual_step_length) <= cutoff1)
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
			= step_length(X_cholesky, dX, dynamical_parameters.solver_parameters.step_length_reduction,
					"run.step.stepLength(XCholesky)", timers);
		dual_step_length
			= step_length(Y_cholesky, dY, dynamical_parameters.solver_parameters.step_length_reduction,
					"run.step.stepLength(YCholesky)", timers);
		if(is_primal_and_dual_feasible)
		{
			primal_step_length = El::Min(primal_step_length, dual_step_length);
			dual_step_length = primal_step_length;
		}
		for(size_t block = 0; block < x.blocks.size(); ++block)
		{
			El::Axpy(primal_step_length, internal_dx.blocks[block], x.blocks[block]);
		}
		// Update the dual point (y, Y) += dualStepLength*(dy, dY)
		for(size_t block = 0; block < dy.blocks.size(); ++block)
		{
			El::Axpy(dual_step_length, internal_dy.blocks[block], y.blocks[block]);
		}
	}
	else
	{
		for(size_t block = 0; block < x.blocks.size(); ++block)
		{
			El::Axpy(primal_step_length, dx.blocks[block], x.blocks[block]);
		}
		// Update the dual point (y, Y) += dualStepLength*(dy, dY)
		for(size_t block = 0; block < dy.blocks.size(); ++block)
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

	/*
	//if (El::Max(primal_step_length, dual_step_length) > 0 && El::Max(primal_step_length, dual_step_length)<0.1)
	//	  mu_direction_mode = 2;
	if (El::Max(primal_step_length, dual_step_length) > 0 && El::Max(primal_step_length, dual_step_length)<0.1) 
		mu_direction_mode = 1;
	else
		mu_direction_mode = 0;
		*/
	//mu_direction_mode = 0;


	step_timer.stop();
} 
