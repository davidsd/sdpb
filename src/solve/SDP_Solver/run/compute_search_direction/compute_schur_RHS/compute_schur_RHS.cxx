#include "../../../../SDP_Solver.hxx"

// Compute the vectors r_x and r_y on the right-hand side of the Schur
// complement equation:
//
// {{S, -B}, {B^T, 0}} . {dx, dy} = {r_x, r_y}
//
// where S = SchurComplement and B = FreeVarMatrix.  Specifically,
//
// r_x[p] = -dualResidues[p] - Tr(A_p Z)              for 0 <= p < P
// r_y[n] = dualObjective[n] - (FreeVarMatrix^T x)_n  for 0 <= n < N
//
// where P = primalObjective.size(), N = dualObjective.size()
//
// Inputs:
// - sdp, an SDP
// - dualResidues, a Vector of length P
// - Z = X^{-1} (PrimalResidues Y - R)
// - x, a vector of length P
// Outputs:
// - r_x, a Vector of length P
// - r_y, a Vector of length N
//

Real bilinear_block_pairing(const Real *v, const int dim, const Matrix &A,
                            const int blockRow, const int blockCol);

void compute_schur_RHS(const SDP &sdp, const Vector &dualResidues,
                       const Block_Diagonal_Matrix &Z, const Vector &x,
                       Vector &r_x, Vector &r_y)
{
  // r_x[p] = -dualResidues[p]
  for(unsigned int p = 0; p < r_x.size(); p++)
    r_x[p] = -dualResidues[p];

  // r_x[p] = -dualResidues[p] - Tr(A_p Z), where A_p are as described
  // in SDP.h.  The trace can be computed in terms of bilinearBases
  // using bilinearBlockPairing.
  for(unsigned int j = 0; j < sdp.dimensions.size(); j++)
    {
      for(auto &t : sdp.constraint_indices[j])
        {
          for(auto &b : sdp.blocks[j])
            {
              const int h = sdp.bilinear_bases[b].rows;
              // Pointer to the k-th column of sdp.bilinearBases[b]
              const Real *q = &sdp.bilinear_bases[b].elements[(t.k) * h];

              r_x[t.p] -= bilinear_block_pairing(q, h, Z.blocks[b], t.r, t.s);
            }
        }
    }

  // r_y[n] = dualObjective[n] - (FreeVarMatrix^T x)_n
  for(unsigned int n = 0; n < sdp.dual_objective_b.size(); n++)
    {
      r_y[n] = sdp.dual_objective_b[n];
      for(int p = 0; p < sdp.free_var_matrix.rows; p++)
        {
          r_y[n] -= sdp.free_var_matrix.elt(p, n) * x[p];
        }
    }
}
