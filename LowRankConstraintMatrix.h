#ifndef SDP_BOOTSTRAP_LOWRANKCONSTRAINTMATRIX_H_
#define SDP_BOOTSTRAP_LOWRANKCONSTRAINTMATRIX_H_

#include <vector>
#include "types.h"

using std::pair;
using std::vector;

// A matrix of the form E_{ij} \otimes v v^T, where E_{ij} is the
// matrix with (i,j)-th entry 1, and the rest zero, and v is a vector.
class RankOneSingleBlockMatrix {
public:
  int i;
  int j;
  vector<Real> vec;

  RankOneSingleBlockMatrix(int i, int j, const vector<Real> &vec):
    i(i), j(j), vec(vec) {}
};

class LowRankConstraintMatrix {
public:
  vector<Real> diagonalPart;
  vector<pair<int, RankOneSingleBlockMatrix> > lowRankBlocks;
  vector<pair<int, DenseMatrix> > denseBlocks;
};

#endif  // SDP_BOOTSTRAP_LOWRANKCONSTRAINTMATRIX_H_
