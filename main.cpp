#include <iostream>
#include <vector>
#include "types.h"
#include "DenseMatrix.h"
#include "BlockDenseMatrix.h"
#include "LowRankConstraintMatrix.h"

using std::cout;
using std::endl;

Real schurComplementPairing(const LowRankConstraintMatrix &F1,
                            const LowRankConstraintMatrix &F2,
                            const BlockDenseMatrix &X,
                            const BlockDenseMatrix &ZInv) {
  Real total = 0;

  for (unsigned int k1 = 0; k1 < F1.lowRankBlocks.size(); k1++) {
    int b1 = F1.lowRankBlocks[k1].first;

    for (unsigned int k2 = 0; k2 < F2.lowRankBlocks.size(); k2++) {
      int b2 = F2.lowRankBlocks[k2].first;

      if (b1 == b2) {
        int i1 = F1.lowRankBlocks[k1].second.i;
        int j1 = F1.lowRankBlocks[k1].second.j;

        int i2 = F2.lowRankBlocks[k2].second.i;
        int j2 = F2.lowRankBlocks[k2].second.j;

        const vector<Real> v = F1.lowRankBlocks[k1].second.vec;
        const vector<Real> u = F2.lowRankBlocks[k2].second.vec;

        total += (X.blocks[b1].blockBilinearForm(j2, i1, u, v) *
                  ZInv.blocks[b1].blockBilinearForm(j1, i2, v, u));

      }
    }
  }
  return total;
        
}

int main(int argc, char** argv) {

  Real x = 1;
  Real y = 2;
  DenseMatrix m(2);
  DenseMatrix n(2);
  m.choleskyDecomposition(n);
  vector<int> sizes;
  sizes.push_back(3);
  BlockDenseMatrix b(vector<Real>(0), sizes);

  cout << x+y << endl;
  for (int i = 0; i < m.dim; i++) {
    for (int j = 0; j < m.dim; j++) {
      cout << b.blocks[0].elem(i,j) << endl;
    }
  }

}
