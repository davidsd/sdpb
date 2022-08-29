#include "../../../sdp_solve/SDP_Solver.hxx"


// Tr(A B), where A and B are symmetric
El::BigFloat frobenius_product_symmetric(const Block_Diagonal_Matrix &A,
	const Block_Diagonal_Matrix &B);

void scale_multiply_add(const El::BigFloat &alpha,
	const Block_Diagonal_Matrix &A,
	const Block_Diagonal_Matrix &B,
	const El::BigFloat &beta, Block_Diagonal_Matrix &C);


// R_error= tr(XY)/X.dim * I - XY

void compute_R_error(const std::size_t &total_psd_rows, const Block_Diagonal_Matrix &X, const Block_Diagonal_Matrix &Y, El::BigFloat & R_error, Timers &timers)
{
  auto &R_error_timer(timers.add_and_start("run.computeRerror"));

  El::BigFloat mu = frobenius_product_symmetric(X, Y) / total_psd_rows;

  Block_Diagonal_Matrix R(X);
  scale_multiply_add(El::BigFloat(-1), X, Y, El::BigFloat(0), R);
  R.add_diagonal(mu);

  R_error = R.max_abs_mpi();

  R_error_timer.stop();
}
