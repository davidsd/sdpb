#ifndef SDP_BOOTSTRAP_BLOCKDENSEMATRIX_H_
#define SDP_BOOTSTRAP_BLOCKDENSEMATRIX_H_

#include <vector>
#include "types.h"
#include "DenseMatrix.h"

using std::vector;

class BlockDenseMatrix {
public:
  vector<Real> diagonalPart;
  vector<DenseMatrix> blocks;

  BlockDenseMatrix(const vector<Real> &diagonalPart, const vector<int> &blockSizes):
    diagonalPart(diagonalPart) {
    for (unsigned int i = 0; i < blockSizes.size(); i++) {
      blocks.push_back(DenseMatrix(blockSizes[i]));
    }
  }

  void inverseCholeskyAndInverse(BlockDenseMatrix &invCholesky,
                                 BlockDenseMatrix &inverse,
                                 BlockDenseMatrix &workSpace) {

    for (unsigned int i = 0; i < diagonalPart.size(); i++) {
      Real d = diagonalPart[i];
      invCholesky.diagonalPart[i] = 1/sqrt(d);
      inverse.diagonalPart[i] = 1/d;
    }

    for (unsigned int b = 0; b < blocks.size(); b++) {
      blocks[b].inverseCholeskyAndInverse(invCholesky.blocks[b], inverse.blocks[b], workSpace.blocks[b]);
    }
  }

};

#endif  // SDP_BOOTSTRAP_BLOCKDENSEMATRIX_H_
