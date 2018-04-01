#include "../../../SDP_Solver.hxx"

// Result = \sum_p a[p] A_p,
//
// where a[p] is a vector of length primalObjective.size() and the
// constraint matrices A_p are given by
//
//   A_(j,r,s,k) = \sum_{b \in blocks[j]}
//                     Block_b(v_{b,k} v_{b,k}^T \otimes E^{rs}),
//
// where v_{b,k} is the k-th column of bilinearBases[b], as described
// in SDP.h.
//
// Inputs: sdp, a
// Output: Result (overwritten)
//

void diagonal_congruence(Real const *d, const Matrix &V, const int blockRow,
                         const int blockCol, Matrix &Result);

void constraint_matrix_weighted_sum(const SDP &sdp, const Vector a,
                                    Block_Diagonal_Matrix &Result)
{
  for(unsigned int j = 0; j < sdp.dimensions.size(); j++)
    {
      const int ej = sdp.degrees[j] + 1;

      // For each j, t points to the first Index_Tuple corresponding to j
      for(auto t(sdp.constraint_indices[j].begin());
          t != sdp.constraint_indices[j].end(); t += ej)
        {
          const int p = t->p;
          const int r = t->r;
          const int s = t->s;
          assert(t->k == 0);

          for(auto &b : sdp.blocks[j])
            {
              // Result.blocks[b]^(r,s) = V diag(a') V^T, where
              // V=sdp.bilinearBases[b], a' denotes the subvector of a
              // corresponding to j, and M^(r,s) denotes the (r,s)-th block
              // of M.
              diagonal_congruence(&a[p], sdp.bilinear_bases[b], r, s,
                                  Result.blocks[b]);

              // Result should be symmetric, so if r != s, we must divide
              // the (r,s)-th block of Result.blocks[b] by 2 and copy its
              // transpose to the (s,r)-th block.
              if(r != s)
                {
                  const int u = sdp.bilinear_bases[b].rows;
                  for(int m = r * u; m < (r + 1) * u; m++)
                    {
                      for(int n = s * u; n < (s + 1) * u; n++)
                        {
                          Result.blocks[b].elt(m, n) /= 2;
                          Result.blocks[b].elt(n, m)
                            = Result.blocks[b].elt(m, n);
                        }
                    }
                }
            }
        }
    }
}
