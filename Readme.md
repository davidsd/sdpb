# SDPB

SDPB is a semidefinite programming solver designed specifically for
application to the conformal bootstrap and featuring arbitrary precision arithmetic.

It solves the following optimization problem

```
minimize \sum_n b_n z_n over vectors z=(z_1,...,z_N) such that

\sum_n z_n M_{k,n}(x)

is a positive semidefinite matrix for all k=1,...,K and all x >= 0.

Here, for each k and n, M_{k,n}(x) is a matrix polynomial in x.  The matrix dimension
of M_{k,n}(x) can depend on k but not n.
```

## Installation and Requirements

SDPB requires

- [Boost C++ Libraries](http://www.boost.org/) (tested with Boost 1.54).

- [The GNU Multiprecision Library](https://gmplib.org/).

To install, you must first edit the `Makefile` to replace `/home/dsd/lib`, `/home/dsd/include` and `/home/dsd/include/boost` with paths pointing to the appropriate places.  Then type `make` to build the `sdpb` executable.

## Usage

Type `sdpb --help` for the syntax and a list of options.

The input format for SDPB is XML-based and described in the manual.  The Mathematica file `Examples.m` includes code to export semidefinite programs in this format, along with some examples.

## Algorithm

SDPB uses a primal-dual interior point method based on...

## The Conformal Bootstrap

## Author

- David Simmons-Duffin (davidsd@gmail.com)

## Attribution

If you use SDPB in work that results in publication, please cite

- D. Simmons-Duffin, "SDPB: A Semidefinite Program Solver for the
  Conformal Bootstrap", arXiv:xxxx.xxxx \[hep-th\].

Depending on how SDPB is used, please also consider the following sources:

The first use of semidefinite programming in the bootstrap:

- D. Poland, D. Simmons-Duffin and A. Vichi, "Carving Out the Space of
  4D CFTs," JHEP 1205, 110 (2012) [arXiv:1109.5176 \[hep-th\]](http://arxiv.org/abs/arXiv:1109.5176).

The generalization of semidefinite programming methods to arbitrary
spacetime dimension:

- F. Kos, D. Poland and D. Simmons-Duffin, "Bootstrapping the O(N)
  Vector Models," JHEP 1406, 091 (2014) [arXiv:1307.6856 \[hep-th\]](http://arxiv.org/abs/arXiv:1307.6856).

The generalization of semidefinite programming methods to arbitrary
systems of correlation functions:

- F. Kos, D. Poland and D. Simmons-Duffin, "Bootstrapping Mixed
  Correlators in the 3D Ising Model," [arXiv:1406.4858 \[hep-th\]](http://arxiv.org/abs/arXiv:1406.4858).

## Acknowledgements

- SDPB Makes extensive use of [MPACK](http://mplapack.sourceforge.net/), the multiple precision linear algebra library written by Nakata Maho.  Several source files from MPACK are included in the SDPB source tree (see the license at the top of those files).

- SDPB uses Lee Thomason's [tinyxml2](http://www.grinninglizard.com/tinyxml2/) library for parsing.

- The design of SDPB was partially based on the solver [SDPA](http://sdpa.sourceforge.net/), without which this software would never have been written.

- Thanks to Filip Kos, David Poland, and Alessandro Vichi for collaboration in developing semidefinite programming methods for the conformal bootstrap.

- Thanks to Amir Ali Ahmadi, Hande Benson, Pablo Parrilo, and Robert Vanderbei for advice and discussions about semidefinite programming.

- Thanks to Luca Iliesiu, Filip Kos, Daliang Li, David Poland, Silviu Pufu, Alessandro Vichi, and Ran Yacoby for assistance testing SDPB.