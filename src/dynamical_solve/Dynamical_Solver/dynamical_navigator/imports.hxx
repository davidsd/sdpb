#include "sdp_solve/sdp_solve.hxx"
#include "dynamical_solve/dynamical_solve.hxx"

//Functions from SDP_Solver used in dynamic_step.cxx
El::BigFloat frobenius_product_symmetric(const Block_Diagonal_Matrix &A,
                                         const Block_Diagonal_Matrix &B);

El::BigFloat predictor_centering_parameter(const Solver_Parameters &parameters,
                                           const bool is_primal_dual_feasible);

El::BigFloat corrector_centering_parameter(
  const Solver_Parameters &parameters, const Block_Diagonal_Matrix &X,
  const Block_Diagonal_Matrix &dX, const Block_Diagonal_Matrix &Y,
  const Block_Diagonal_Matrix &dY, const El::BigFloat &mu,
  const bool is_primal_dual_feasible, const size_t &total_num_rows);

void cholesky_decomposition(const Block_Diagonal_Matrix &A,
                            Block_Diagonal_Matrix &L);

void cholesky_solve(const Block_Diagonal_Matrix &ACholesky,
                    Block_Diagonal_Matrix &X);

void constraint_matrix_weighted_sum(const Block_Info &block_info,
                                    const SDP &sdp, const Block_Vector &a,
                                    Block_Diagonal_Matrix &result);

El::BigFloat
step_length(const Block_Diagonal_Matrix &MCholesky,
            const Block_Diagonal_Matrix &dM, const El::BigFloat &gamma,
            const std::string &timer_name, Timers &timers);

void initialize_schur_complement_solver(
  const Block_Info &block_info, const SDP &sdp,
  const std::array<
    std::vector<std::vector<std::vector<El::DistMatrix<El::BigFloat>>>>, 2>
    &A_X_inv,
  const std::array<
    std::vector<std::vector<std::vector<El::DistMatrix<El::BigFloat>>>>, 2>
    &A_Y,
  const El::Grid &group_grid, Block_Diagonal_Matrix &schur_complement_cholesky,
  Block_Matrix &schur_off_diagonal,
  BigInt_Shared_Memory_Syrk_Context &bigint_syrk_context,
  El::DistMatrix<El::BigFloat> &Q, Timers &timers,
  El::Matrix<int32_t> &block_timings_ms);

void Axpy(const El::BigFloat &alpha, const SDP &new_sdp, SDP &delta_sdp);

// C := alpha*A*B + beta*C
void scale_multiply_add(const El::BigFloat &alpha,
                        const Block_Diagonal_Matrix &A,
                        const Block_Diagonal_Matrix &B,
                        const El::BigFloat &beta, Block_Diagonal_Matrix &C);

// C := A B
inline void multiply(const Block_Diagonal_Matrix &A,
                     const Block_Diagonal_Matrix &B, Block_Diagonal_Matrix &C)
{
  scale_multiply_add(El::BigFloat(1), A, B, El::BigFloat(0), C);
}
