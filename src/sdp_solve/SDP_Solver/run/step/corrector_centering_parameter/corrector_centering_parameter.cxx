#include "sdp_solve/SDP_Solver.hxx"

// (X + dX) . (Y + dY), where X, dX, Y, dY are symmetric
// BlockDiagonalMatrices and '.' is the Frobenius product.
//
El::BigFloat frobenius_product_of_sums(const Block_Diagonal_Matrix &X,
                                       const Block_Diagonal_Matrix &dX,
                                       const Block_Diagonal_Matrix &Y,
                                       const Block_Diagonal_Matrix &dY);

// Centering parameter \beta_c for the corrector step
El::BigFloat corrector_centering_parameter(
  const Solver_Parameters &parameters, const Block_Diagonal_Matrix &X,
  const Block_Diagonal_Matrix &dX, const Block_Diagonal_Matrix &Y,
  const Block_Diagonal_Matrix &dY, const El::BigFloat &mu,
  const bool is_primal_dual_feasible, const size_t &total_psd_rows)
{
  El::BigFloat r
    = frobenius_product_of_sums(X, dX, Y, dY) / (mu * total_psd_rows);
  El::BigFloat beta = r < 1 ? r * r : r;

  if(is_primal_dual_feasible)
    {
      return Min(Max(parameters.feasible_centering_parameter, beta),
                 El::BigFloat(1));
    }
  else
    {
      return Max(parameters.infeasible_centering_parameter, beta);
    }
}
