#include "../../SDP_Solver.hxx"

// dualResidues[p] = primalObjective[p] - Tr(A_p Y) - (FreeVarMatrix y)_p,
// for 0 <= p < primalObjective.size()
//
// The pairings Tr(A_p Y) can be written in terms of BilinearPairingsY:
//
//   Tr(A_(j,r,s,k) Y) = \sum_{b \in blocks[j]}
//                       (1/2) (BilinearPairingsY_{ej r + k, ej s + k} +
//                              swap (r <-> s))
// where ej = d_j + 1.
//
// Inputs: sdp, y, BilinearPairingsY
// Output: dualResidues (overwriten)
//
void compute_dual_residues(const SDP &sdp, const Vector &y,
                           const Block_Diagonal_Matrix &bilinear_pairings_Y,
                           Vector &dual_residues)
{
  for(unsigned int j = 0; j < sdp.dimensions.size(); j++)
    {
      const int ej = sdp.degrees[j] + 1;

      for(auto &t : sdp.constraint_indices[j])
        {
          const int p = t.p;
          const int ej_r = t.r * ej;
          const int ej_s = t.s * ej;
          const int k = t.k;

          // dualResidues[p] = -Tr(A_p Y)
          dual_residues[p] = 0;
          for(auto &b : sdp.blocks[j])
            {
              dual_residues[p]
                -= bilinear_pairings_Y.blocks[b].elt(ej_r + k, ej_s + k);
              dual_residues[p]
                -= bilinear_pairings_Y.blocks[b].elt(ej_s + k, ej_r + k);
            }
          dual_residues[p] /= 2;

          // dualResidues[p] = -Tr(A_p Y) - (FreeVarMatrix y)_p
          for(int n = 0; n < sdp.free_var_matrix.cols; n++)
            dual_residues[p] -= sdp.free_var_matrix.elt(p, n) * y[n];

          // dualResidues[p] = primalObjective[p] - Tr(A_p Y) - (FreeVarMatrix
          // y)_p
          dual_residues[p] += sdp.primal_objective_c[p];
        }
    }
}
