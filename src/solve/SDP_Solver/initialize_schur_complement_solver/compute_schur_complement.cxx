#include "../../SDP_Solver.hxx"

// Compute the SchurComplement matrix using BilinearPairingsXInv and
// BilinearPairingsY and the formula
//
//   S_{(j,r1,s1,k1), (j,r2,s2,k2)} = \sum_{b \in blocks[j]}
//          (1/4) (BilinearPairingsXInv_{ej s1 + k1, ej r2 + k2}*
//                 BilinearPairingsY_{ej s2 + k2, ej r1 + k1} +
//                 swaps (r1 <-> s1) and (r2 <-> s2))
//
// where ej = d_j + 1.
//
// Inputs: sdp, BilinearPairingsXInv, BilinearPairingsY
// Output: SchurComplement (overwritten)
//
void compute_schur_complement(
  const SDP &sdp, const Block_Diagonal_Matrix &bilinear_pairings_X_inv,
  const Block_Diagonal_Matrix &bilinear_pairings_Y,
  Block_Diagonal_Matrix &schur_complement)
{
  for(unsigned int j = 0; j < sdp.dimensions.size(); j++)
    {
      const int ej = sdp.degrees[j] + 1;

      for(unsigned int u1 = 0; u1 < sdp.constraint_indices[j].size(); u1++)
        {
          const int ej_r1 = sdp.constraint_indices[j][u1].r * ej;
          const int ej_s1 = sdp.constraint_indices[j][u1].s * ej;
          const int k1 = sdp.constraint_indices[j][u1].k;

          for(unsigned int u2 = 0; u2 <= u1; u2++)
            {
              const int ej_r2 = sdp.constraint_indices[j][u2].r * ej;
              const int ej_s2 = sdp.constraint_indices[j][u2].s * ej;
              const int k2 = sdp.constraint_indices[j][u2].k;

              Real tmp = 0;
              for(auto &b : sdp.blocks[j])
                {
                  tmp += (bilinear_pairings_X_inv.blocks[b].elt(ej_s1 + k1,
                                                                ej_r2 + k2)
                            * bilinear_pairings_Y.blocks[b].elt(ej_s2 + k2,
                                                                ej_r1 + k1)
                          + bilinear_pairings_X_inv.blocks[b].elt(ej_r1 + k1,
                                                                  ej_r2 + k2)
                              * bilinear_pairings_Y.blocks[b].elt(ej_s2 + k2,
                                                                  ej_s1 + k1)
                          + bilinear_pairings_X_inv.blocks[b].elt(ej_s1 + k1,
                                                                  ej_s2 + k2)
                              * bilinear_pairings_Y.blocks[b].elt(ej_r2 + k2,
                                                                  ej_r1 + k1)
                          + bilinear_pairings_X_inv.blocks[b].elt(ej_r1 + k1,
                                                                  ej_s2 + k2)
                              * bilinear_pairings_Y.blocks[b].elt(ej_r2 + k2,
                                                                  ej_s1 + k1))
                         / 4;
                }
              schur_complement.blocks[j].elt(u1, u2) = tmp;
              if(u2 != u1)
                {
                  schur_complement.blocks[j].elt(u2, u1) = tmp;
                }
            }
        }
    }
}
