#include "../../../sdp_solve.hxx"
#include "../../../dynamical_solve.hxx"

//Compute dx and dy of the central sdp as the standard sdp_solver does. 
//Correspond to - H^-1_xx Del_p L_mu in Eq(13).
//Return: void.
//Update dx, dy.

void scale_multiply_add(const El::BigFloat &alpha,
                        const Block_Diagonal_Matrix &A,
                        const Block_Diagonal_Matrix &B,
                        const El::BigFloat &beta, Block_Diagonal_Matrix &C);


inline void multiply(const Block_Diagonal_Matrix &A,
                     const Block_Diagonal_Matrix &B, Block_Diagonal_Matrix &C)
{
  scale_multiply_add(El::BigFloat(1), A, B, El::BigFloat(0), C);
}

void cholesky_solve(const Block_Diagonal_Matrix &ACholesky,
                    Block_Diagonal_Matrix &X);


void compute_schur_RHS(const Block_Info &block_info, const SDP &sdp,
                       const Block_Vector &dual_residues,
                       const Block_Diagonal_Matrix &Z,
                       Block_Vector &dx);

void solve_schur_complement_equation(
  const Block_Diagonal_Matrix &schur_complement_cholesky,
  const Block_Matrix &schur_off_diagonal,
  const El::DistMatrix<El::BigFloat> &Q, Block_Vector &dx, Block_Vector &dy);

void internal_predictor_direction(
  const Block_Info &block_info, const SDP &sdp, const Dynamical_Solver &solver,
  const Block_Diagonal_Matrix &schur_complement_cholesky,
  const Block_Matrix &schur_off_diagonal,
  const Block_Diagonal_Matrix &X_cholesky, const El::BigFloat beta,
  const El::BigFloat &mu, const Block_Vector &primal_residue_p,
  const El::DistMatrix<El::BigFloat> &Q, Block_Vector &grad_x, Block_Vector &grad_y, Block_Vector &dx, Block_Vector &dy, Block_Diagonal_Matrix &R)
{

//  Block_Diagonal_Matrix R(solver.X);

  scale_multiply_add(El::BigFloat(-1), solver.X, solver.Y, El::BigFloat(0), R);
  R.add_diagonal(beta * mu);

  Block_Diagonal_Matrix Z(solver.X);
  multiply(solver.primal_residues, solver.Y, Z);
  Z -= R;
  cholesky_solve(X_cholesky, Z);
  Z.symmetrize();

  compute_schur_RHS(block_info, sdp, solver.dual_residues, Z, dx);
  dy=primal_residue_p;

  grad_x = dx;
  grad_y = dy; 

  solve_schur_complement_equation(schur_complement_cholesky,
                                  schur_off_diagonal, Q, dx, dy);
}  


