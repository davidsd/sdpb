#include <vector>
#include <boost/filesystem.hpp>
#include <boost/filesystem/fstream.hpp>
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


void compute_schur_RHS(const Block_Info &block_info, const SDP &sdp,
	const Block_Vector &dual_residues,
	const Block_Diagonal_Matrix &Z,
	Block_Vector &dx);

void solve_schur_complement_equation(
	const Block_Diagonal_Matrix &schur_complement_cholesky,
	const Block_Matrix &schur_off_diagonal,
	const El::DistMatrix<El::BigFloat> &Q, Block_Vector &dx, Block_Vector &dy);


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

void compute_R_error(const std::size_t &total_psd_rows, const Block_Diagonal_Matrix &X, const Block_Diagonal_Matrix &Y,
	Block_Diagonal_Matrix &R, El::BigFloat & R_error, El::BigFloat & mu, Timers &timers);

void compute_R_error(const std::size_t &total_psd_rows,
	const Block_Diagonal_Matrix &X, const Block_Diagonal_Matrix &Y, El::BigFloat & R_error, Timers &timers);


void compute_corrector_errors(
	const Block_Info &block_info, const SDP &sdp,
	const std::size_t &total_psd_rows,

	const Block_Vector &x_const, const Block_Vector &dx_const,
	const Block_Vector &y_const, const Block_Vector &dy_const,
	const Block_Diagonal_Matrix &X_const, const Block_Diagonal_Matrix &dX_const,
	const Block_Diagonal_Matrix &Y_const, const Block_Diagonal_Matrix &dY_const,

	const El::BigFloat &primal_step_length, const El::BigFloat &dual_step_length,

	El::BigFloat &primal_error_P,
	El::BigFloat &primal_error_p,
	El::BigFloat &dual_error,
	El::BigFloat &R_error,
	El::BigFloat &mu,

	Timers &timers);

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

void compute_corrector_R(
	const Block_Info &block_info, const std::size_t &total_psd_rows,

	const Block_Vector &x_const, const Block_Vector &dx_const,
	const Block_Vector &y_const, const Block_Vector &dy_const,
	const Block_Diagonal_Matrix &X_const, const Block_Diagonal_Matrix &dX_const,
	const Block_Diagonal_Matrix &Y_const, const Block_Diagonal_Matrix &dY_const,

	El::BigFloat &R_error,
	El::BigFloat &mu,

	Timers &timers);


////////////////////////////////////
// to Aike : maybe we should package those functions as class member functions?
// what's the difference between El::Zero and El::Zeros
void Block_Vector_add(const Block_Vector &dv, Block_Vector &v)
{
	for (size_t block = 0; block < dv.blocks.size(); ++block)
		El::Axpy(1, dv.blocks[block], v.blocks[block]);
}

void Block_Vector_zero(Block_Vector &v)
{
	for (size_t block = 0; block < v.blocks.size(); ++block)
		El::Zero(v.blocks[block]);
}

void Block_Vector_negative(Block_Vector &v)
{
	for (size_t block = 0; block < v.blocks.size(); ++block)
		v.blocks[block] *= -1;
}


void Block_Diagonal_Matrix_zero(Block_Diagonal_Matrix &C)
{
	for (size_t block = 0; block < C.blocks.size(); ++block)
		El::Zero(C.blocks[block]);
}

void Block_Diagonal_Matrix_negative(Block_Diagonal_Matrix &C)
{
	for (size_t block = 0; block < C.blocks.size(); ++block)
		C.blocks[block] *= -1;
}
////////////////////////////////////



void Dynamical_Solver::external_corrector_step(const Dynamical_Solver_Parameters &dynamical_parameters,
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
	Timers &timers)
{
	Block_Diagonal_Matrix primal_residues_0(X);
	Block_Vector dual_residues_0(block_info.schur_block_sizes(), block_info.block_indices,
		block_info.num_points.size(), grid);
	Block_Vector primal_residue_p_0(y);
	Block_Diagonal_Matrix R_0(X);
	El::BigFloat primal_error_P_0, primal_error_p_0, dual_error_0, R_error_0;
	El::BigFloat coit_mu;

	// compute P_0, p_0, d_0, R_0
	compute_dual_residues_and_error(block_info, old_sdp, y, A_Y, dual_residues_0,
		dual_error_0, timers);
	compute_primal_residues_and_error_P_Ax_X(
		block_info, old_sdp, x, X, primal_residues_0, primal_error_P_0, timers);
	compute_primal_residues_and_error_p_b_Bx(
		block_info, old_sdp, x, primal_residue_p_0, primal_error_p_0);
	compute_corrector_R(block_info, total_psd_rows,
		x, dx, y, dy, X, dX, Y, dY,
		1, 1, R_0, R_error_0, coit_mu, timers);

	scale_multiply_add(El::BigFloat(-1), X, Y, El::BigFloat(0), R_0);
	R_0.add_diagonal(coit_mu);

	Block_Diagonal_Matrix R_old(X); El::BigFloat R_err_old, coit_mu_old;
	compute_R_error(total_psd_rows, X, Y, R_old, R_err_old, coit_mu_old, timers);
	
	if (El::mpi::Rank() == 0)std::cout << "before external correcter step, R_err(X)=" << R_err_old
		<< " mu=" << coit_mu_old
		<< " R_err(X+dX)=" << R_error_0
		<< " mu=" << coit_mu << "\n" << std::flush;
		

	SDP d_sdp(new_sdp);
	Axpy(El::BigFloat(-1), old_sdp, d_sdp);

	// compute P, p, d, R
	Block_Diagonal_Matrix primal_residues(primal_residues_0);
	Block_Vector dual_residues(dual_residues_0);
	Block_Vector primal_residue_p(primal_residue_p_0);
	Block_Diagonal_Matrix R(R_0);

	compute_corrector_residue_shift(block_info,
		primal_residue_p_0, dual_residues_0, R_0,
		primal_residue_p, dual_residues, R,
		dx, dy, dX, dY,
		d_sdp);

	// compute RHS
	Block_Diagonal_Matrix Z(X);
	multiply(primal_residues, Y, Z);
	Z -= R;
	cholesky_solve(X_cholesky, Z);
	Z.symmetrize();

	compute_schur_RHS(block_info, old_sdp, dual_residues, Z, dx);
	dy = primal_residue_p;

	solve_schur_complement_equation(schur_complement_cholesky,
		schur_off_diagonal, Q, dx, dy);

	compute_dXdY(true, block_info, old_sdp, grid,
		X_cholesky, Y_cholesky, timers,
		dx, dy, dX, dY, R,
		primal_step_length, dual_step_length, step_length_reduction);

	///////// test ///////////////
	El::BigFloat error_P, error_p, error_d, error_R;
	compute_corrector_errors(block_info, new_sdp, total_psd_rows,
		x, dx, y, dy, X, dX, Y, dY,
		1, 1,
		error_P, error_p, error_d, error_R, coit_mu,
		timers);
	 
	if (El::mpi::Rank() == 0)std::cout << "external corrector : "
		<< "R=" << error_R 
		<< " mu=" << coit_mu
		<< " step-len : ("
		<< primal_step_length << ", " << dual_step_length << ")\n" << std::flush;

	current_error_R = error_R;
	current_mu = coit_mu;

	return;
}

void read_Schur_data(boost::filesystem::path schur_path, const Block_Info &block_info,
	Block_Diagonal_Matrix &schur_complement_cholesky,
	Block_Matrix &schur_off_diagonal,
	El::DistMatrix<El::BigFloat> &Q)
{
	if (boost::filesystem::exists(schur_path / "Q_cholesky.txt"))
	{
		for (size_t block = 0; block != block_info.block_indices.size(); ++block)
		{
			size_t block_index(block_info.block_indices.at(block));
			read_text_block(schur_complement_cholesky.blocks.at(block),
				schur_path, "schur_complement_cholesky_",
				block_index);
			read_text_block(schur_off_diagonal.blocks.at(block),
				schur_path, "schur_off_diagonal_", block_index);
		}
		read_text_block(Q, schur_path / "Q_cholesky.txt");
	}
}


bool load_binary_checkpoint(const boost::filesystem::path &checkpoint_directory,
	Block_Vector & dx, Block_Vector & dy,
	Block_Diagonal_Matrix & dX, Block_Diagonal_Matrix & dY,
	const Verbosity &verbosity);

extern Verbosity verbosity_save;

void Dynamical_Solver::external_corrector_run(const Dynamical_Solver_Parameters &dynamical_parameters,
	const SDP &new_sdp,
	const Block_Info &block_info, const El::Grid &grid,
	Block_Vector & dx, Block_Vector & dy, Block_Diagonal_Matrix & dX, Block_Diagonal_Matrix & dY,
	El::BigFloat &primal_step_length, El::BigFloat &dual_step_length, El::BigFloat &step_length_reduction,
	Timers &timers)
{
	auto psd_sizes(block_info.psd_matrix_block_sizes());
	std::size_t total_psd_rows(
		std::accumulate(psd_sizes.begin(), psd_sizes.end(), size_t(0)));

	SDP old_sdp(dynamical_parameters.old_sdp_path, block_info, grid);

	Block_Diagonal_Matrix schur_complement_cholesky(
		block_info.schur_block_sizes(), block_info.block_indices,
		block_info.num_points.size(), grid);
	Block_Matrix schur_off_diagonal(new_sdp.free_var_matrix);
	El::DistMatrix<El::BigFloat> Q(new_sdp.dual_objective_b.Height(),
		new_sdp.dual_objective_b.Height());

	Block_Diagonal_Matrix X_cholesky(X), Y_cholesky(X);
	std::array<
		std::vector<std::vector<std::vector<El::DistMatrix<El::BigFloat>>>>, 2>
		A_X_inv;
	std::array<
		std::vector<std::vector<std::vector<El::DistMatrix<El::BigFloat>>>>, 2>
		A_Y;

	cholesky_decomposition(X, X_cholesky);
	cholesky_decomposition(Y, Y_cholesky);

	compute_bilinear_pairings(block_info, X_cholesky, Y, old_sdp.bases_blocks,
		A_X_inv, A_Y, timers);

	read_Schur_data(dynamical_parameters.old_schur_path / (".schur"), block_info,
		schur_complement_cholesky, schur_off_diagonal, Q);

	//Block_Vector_zero(dx);Block_Vector_zero(dy);Block_Diagonal_Matrix_zero(dX);Block_Diagonal_Matrix_zero(dY);

	load_binary_checkpoint(dynamical_parameters.old_schur_path / (".ext_ck"),
		dx, dy, dX, dY, verbosity_save);

	Block_Vector dx_last(dx);
	Block_Vector dy_last(dy);
	Block_Diagonal_Matrix dX_last(dX);
	Block_Diagonal_Matrix dY_last(dY);
	El::BigFloat error_R_last, mu_last, current_error_R, current_mu, primal_step_length_last, dual_step_length_last;

	dx_last = dx;
	dy_last = dy;
	dX_last = dX;
	dY_last = dY;
	compute_corrector_R(block_info, total_psd_rows,
		x, dx, y, dy, X, dX, Y, dY,
		error_R_last, mu_last, timers);

	for (int i = 0; i < 5; i++)
	{
		external_corrector_step(dynamical_parameters, old_sdp, new_sdp,
			block_info, grid, total_psd_rows,
			X_cholesky, Y_cholesky, A_X_inv, A_Y,
			dx, dy,
			dX, dY,
			schur_complement_cholesky, schur_off_diagonal, Q,
			current_error_R, current_mu,
			primal_step_length, dual_step_length, step_length_reduction,
			timers);

		if (El::mpi::Rank() == 0) std::cout
			<< "external corrector : R_error=" << current_error_R
			<< " mu=" << current_mu << "\n" << std::flush;

		if (current_error_R > error_R_last)break;

		dx_last = dx;
		dy_last = dy;
		dX_last = dX;
		dY_last = dY;
		error_R_last = current_error_R;
	}

	primal_step_length = step_length(X_cholesky, dX_last, El::BigFloat(0.9), "run.step.stepLength(XCholesky)", timers);
	dual_step_length = step_length(Y_cholesky, dY_last, El::BigFloat(0.9), "run.step.stepLength(YCholesky)", timers);

	execute_step(dx_last, dy_last, dX_last, dY_last, primal_step_length, dual_step_length);
}




void compute_By(const Block_Info &block_info, const SDP &sdp, const Block_Vector &y, Block_Vector &By)
{
	auto By_block(By.blocks.begin());
	auto primal_objective_c_block(sdp.primal_objective_c.blocks.begin());
	auto y_block(y.blocks.begin());
	auto free_var_matrix_block(sdp.free_var_matrix.blocks.begin());

	for (auto &block_index : block_info.block_indices)
	{
		// By = 0
		Zero(*By_block);

		// By += FreeVarMatrix * y
		Gemm(El::Orientation::NORMAL, El::Orientation::NORMAL, El::BigFloat(1),
			*free_var_matrix_block, *y_block, El::BigFloat(1),
			*By_block);

		++y_block;
		++free_var_matrix_block;
		++By_block;
	}
}


void compute_db_minus_dBnewx(const Block_Info &block_info,
	const SDP &d_sdp,
	const Block_Vector &x0, const Block_Vector &dx,
	Block_Vector &rslt)
{
	Block_Vector x(x0);
	Block_Vector_add(dx, x);

	auto free_var_matrix_block(d_sdp.free_var_matrix.blocks.begin());
	auto x_block(x.blocks.begin());
	auto rslt_block(rslt.blocks.begin());

	for (auto &block_index : block_info.block_indices)
	{
		El::Gemv(El::OrientationNS::TRANSPOSE, El::BigFloat(-1),
			*free_var_matrix_block, *x_block, El::BigFloat(0),
			*rslt_block);

		// The total primal error is the sum of all of the different
		// blocks.  So to prevent double counting, only add
		// dual_objective_b to one of the residue blocks.
		if (block_index == 0)
		{
			El::Axpy(El::BigFloat(1), d_sdp.dual_objective_b, *rslt_block);
		}

		++free_var_matrix_block;
		++x_block;
		++rslt_block;
	}
} 

void scale_multiply_add(const El::BigFloat &alpha,
	const Block_Diagonal_Matrix &A,
	const Block_Diagonal_Matrix &B,
	const El::BigFloat &beta, Block_Diagonal_Matrix &C);



void Dynamical_Solver::compute_corrector_residue_shift(const Block_Info &block_info,
	Block_Vector &primal_residue_p_0, Block_Vector &dual_residues_0, Block_Diagonal_Matrix &R_0,
	Block_Vector &primal_residue_p, Block_Vector &dual_residues, Block_Diagonal_Matrix &R,
	Block_Vector & dx, Block_Vector & dy, Block_Diagonal_Matrix & dX, Block_Diagonal_Matrix & dY,
	const SDP &d_sdp)
{
	Block_Vector dBy_plus_dBdy(x);
	Block_Vector newy(y);
	Block_Vector_add(dy, newy);
	Block_Vector_negative(newy);
	// newy = -y-dy

	// d=d0+dc-dB.y-dB.dy
	compute_By(block_info, d_sdp, newy, dBy_plus_dBdy);
	dual_residues = dual_residues_0;
	Block_Vector_add(d_sdp.primal_objective_c, dual_residues);
	Block_Vector_add(dBy_plus_dBdy, dual_residues);

	Block_Vector p_shift(y);
	compute_db_minus_dBnewx(block_info, d_sdp, x, dx, p_shift);
	primal_residue_p = primal_residue_p_0;
	Block_Vector_add(p_shift, primal_residue_p);

	Block_Diagonal_Matrix R_shift(dX);
	scale_multiply_add(El::BigFloat(-1), dX, dY, El::BigFloat(0), R_shift);
	R = R_0;
	R += R_shift;

	return;
}
