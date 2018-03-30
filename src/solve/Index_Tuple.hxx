#pragma once

// Instead of working with the 1 to 1 correspondence p <-> (j,r,s,k),
// it is convenient to collect tuples with the same index j. This
// gives a 1-to-many correspondence
//
// j <-> { list of IndexTuple(p,r,s,k) }
//
// IndexTuple is simply a named version of the 4-tuple (p,r,s,k).  A
// given IndexTuple uniquely specifies a constraint matrix A_p.
//
class Index_Tuple
{
public:
  int p; // overall index of the constraint
  int r; // first index for E^{rs}
  int s; // second index for E^{rs}
  int k; // index for v_{b,k}
  Index_Tuple(int p, int r, int s, int k) : p(p), r(r), s(s), k(k) {}
  Index_Tuple() {}
};
