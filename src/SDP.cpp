#include "SDP.h"

SDP bootstrapSDP(const Vector &objective,
                 const Real &objectiveConst,
                 const vector<SampledMatrixPolynomial> &sampledMatrixPols) {
  SDP sdp;
  sdp.dualObjective  = objective;
  sdp.objectiveConst = objectiveConst;

  for (vector<SampledMatrixPolynomial>::const_iterator s = sampledMatrixPols.begin();
       s != sampledMatrixPols.end();
       s++) {
    sdp.dimensions.push_back(s->dim);
    sdp.degrees.push_back(s->degree);
    sdp.primalObjective.insert(sdp.primalObjective.end(),
                               s->constraintConstants.begin(),
                               s->constraintConstants.end());
  }
  sdp.FreeVarMatrix = Matrix(sdp.primalObjective.size(), sdp.dualObjective.size());

  int p = 0;
  for (vector<SampledMatrixPolynomial>::const_iterator s = sampledMatrixPols.begin();
       s != sampledMatrixPols.end();
       s++) {

    vector<int> blocks;
    for (vector<Matrix>::const_iterator b = s->bilinearBases.begin();
         b != s->bilinearBases.end();
         b++) {
      assert(b->cols == s->degree + 1);
      blocks.push_back(sdp.bilinearBases.size());
      sdp.bilinearBases.push_back(*b);
    }
    sdp.blocks.push_back(blocks);

    for (int k = 0; k < s->constraintMatrix.rows; k++, p++)
      for (int n = 0; n < s->constraintMatrix.cols; n++)
        sdp.FreeVarMatrix.elt(p, n) = s->constraintMatrix.elt(k, n);
  }
  assert(p == (int)sdp.primalObjective.size());

  sdp.initializeConstraintIndices();
  return sdp;
}
