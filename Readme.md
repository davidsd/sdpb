## Contents

* [SDPB](#sdpb)
* [Installation and Requirements](#installation-and-requirements)
  * [Linux Installation](#linux-installation)
  * [Mac OS X Installation](#mac-os-x-installation)
  * [Windows Installation](#windows-installation)
* [Usage](#usage)
* [Changelog](#changelog)
* [Attribution](#attribution)
* [Acknowledgements](#acknowledgements)
* [Works Using SDPB](#works-using-sdpb)

# SDPB

SDPB is an open-source, arbitrary-precision, parallelized semidefinite
program solver, designed for the conformal bootstrap.

It solves the following type of optimization problem

```
maximize:  b_0 + \sum_n b_n y_n over (y_1,...,y_N),

such that: M_{0j}(x) + \sum_n y_n M_{nj}(x) is positive semidefinite
           for all x >= 0 and 1 <= j <= J,

where each M_{nj}(x) is a polynomial matrix in x.
```

For more information, see [A Semidefinite Program Solver for the Conformal Bootstrap](http://arxiv.org/abs/1502.02033)
and [the manual](https://github.com/davidsd/sdpb/blob/master/docs/SDPB-Manual.pdf).

Author: David Simmons-Duffin (davidsd@gmail.com). As of February 2015, I am
supported by DOE grant number DE-SC0009988 and a William D. Loughlin Membership
at the [Institute for Advanced Study](http://sns.ias.edu).

## Installation and Requirements

SDPB requires

- [Boost C++ Libraries](http://www.boost.org/) (version 1.54 or later).

- [The GNU Multiprecision Library](https://gmplib.org/).

System-specific installation instructions are below.  Please use the
issue tracker if you encounter installation problems. For those with
experience packaging software, I'd appreciate help making SDPB easier
to install.

### Linux Installation

SDPB has been tested on Red Hat Linux. (Many thanks to Chris Beem for helping with installation instructions.) To install,

1. Download Boost and GMP from the links above. Install GMP with the option `--enable-cxx` added to `./configure`. Install Boost.

2. Edit the `Makefile` to define the variables `GMPINCLUDEDIR`,
`BOOSTINCLUDEDIR`, and `LIBDIR.` Ensure `LIBDIR` is in your `LD_LIBRARY_PATH`. 

3. Type `make` to build the `sdpb` executable.

### Mac OS X Installation

The following instructions have been tested on Mac OS 10.10 Yosemite.  (Many thanks to Ying Lin.)

1. Install Homebrew and `gcc-4.9` (or later), for instance by running the following commands

        # Install Homebrew
        ruby -e "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/master/install)"

        # Update Ruby
        brew install ruby

        # Install the latest gcc as well as its dependencies
        # the option --without-multilib avoids a bug in OpenMP support
        brew install gcc --without-multilib

2. Make `/usr/local/bin/g++-4.9` (or whatever version you have) the default compiler by renaming `gcc` and `g++` in `/usr/bin` and creating symlinks

        ln -s /usr/local/bin/g++-4.9 /usr/local/bin/g++
        ln -s /usr/local/bin/gcc-4.9 /usr/local/bin/gcc

3. Unfortunately, Homebrew's versions of GMP and Boost will not work -- they must be compiled from source. Download the latest GMP from [the GMP website](https://gmplib.org/). Upack the tarball (you may need `lzip` which you can install with `brew install lzip`) and `cd` to the `gmp` directory.  Run

        ./configure --enable-cxx
        make
        make install

4. Download Boost from [the Boost website](http://www.boost.org/).  Unpack the tarball and `cd` to the `boost` directory. Run

        ./bootstrap.sh
        ./b2
        sudo ./b2 install
        
   (Note that `bootstrap.sh` above is just an installation script and has absolutely nothing
   to do with the conformal bootstrap -- lots of people like the name "bootstrap"!)
        
5. Type `make` in the `sdpb` directory to compile the `sdpb` executable. If you installed any of the
above software in custom locations, you'll have to modify variables in the
`Makefile` as described in the Linux instructions.

### Windows Installation

Details on Windows installation using Cygwin can be found [here](Readme_Win.md).

## Usage

Type `sdpb --help` for the syntax and a list of options.

The input format for SDPB is XML-based and described in the manual.
The Mathematica file `SDPB.m` includes code to export semidefinite
programs in this format, along with some examples.

## Changelog

- 5/29/15: Fixed a bug in `SDP.cpp` that caused SDPB to incorrectly load
  matrices with dimensions larger than 2x2. (Thanks to Filip Kos for
  the fix.)
- 10/20/15: Fixed a Windows incompatibility in the parser ([Pull request](https://github.com/davidsd/sdpb/pull/8))

## Attribution

If you use SDPB in work that results in publication, please cite

- D. Simmons-Duffin, "A Semidefinite Program Solver for the
  Conformal Bootstrap", [arXiv:1502.02033 \[hep-th\]](http://arxiv.org/abs/1502.02033).

Depending on how SDPB is used, please also consider the following sources:

The first use of semidefinite programming in the bootstrap:

- D. Poland, D. Simmons-Duffin and A. Vichi, "Carving Out the Space of
  4D CFTs," JHEP 1205, 110 (2012) [arXiv:1109.5176 \[hep-th\]](http://arxiv.org/abs/1109.5176).

The generalization of semidefinite programming methods to arbitrary
spacetime dimension:

- F. Kos, D. Poland and D. Simmons-Duffin, "Bootstrapping the O(N)
  Vector Models," JHEP 1406, 091 (2014) [arXiv:1307.6856 \[hep-th\]](http://arxiv.org/abs/1307.6856).

The generalization of semidefinite programming methods to arbitrary
systems of correlation functions:

- F. Kos, D. Poland and D. Simmons-Duffin, "Bootstrapping Mixed
  Correlators in the 3D Ising Model," [arXiv:1406.4858 \[hep-th\]](http://arxiv.org/abs/1406.4858).

## Acknowledgements

- SDPB Makes extensive use of [MPACK](http://mplapack.sourceforge.net/), the multiple precision linear algebra library written by Nakata Maho.  Several source files from MPACK are included in the SDPB source tree (see the license at the top of those files).

- SDPB uses Lee Thomason's [tinyxml2](http://www.grinninglizard.com/tinyxml2/) library (included in the source tree) for parsing.

- The design of SDPB was partially based on the solvers [SDPA](http://sdpa.sourceforge.net/) and SDPA-GMP, which were essential sources of inspiration and examples.

- Thanks to Filip Kos, David Poland, and Alessandro Vichi for collaboration in developing semidefinite programming methods for the conformal bootstrap and assistance testing SDPB.

- Thanks to Amir Ali Ahmadi, Hande Benson, Pablo Parrilo, and Robert Vanderbei for advice and discussions about semidefinite programming.

- Thanks also to Noah Stein, who first suggested the idea of semidefinite programming to me in [this Math Overflow question](http://mathoverflow.net/questions/33242/continuous-linear-programming-estimating-a-solution).

## Works Using SDPB

- F. Kos, D. Poland, D. Simmons-Duffin and A. Vichi,
  "Bootstrapping the O(N) Archipelago,"
  [arXiv:1504.07997 [hep-th]](http://arxiv.org/abs/1504.07997).

- S. M. Chester, S. Giombi, L. V. Iliesiu, I. R. Klebanov, S. S. Pufu and R. Yacoby,
  "Accidental Symmetries and the Conformal Bootstrap,"
  [arXiv:1507.04424 [hep-th]](http://arxiv.org/abs/1507.04424).

- C. Beem, M. Lemos, L. Rastelli and B. C. van Rees,
  "The (2,0) superconformal bootstrap,"
  [arXiv:1507.05637 [hep-th]](http://arxiv.org/abs/1507.05637).

- L. Iliesiu, F. Kos, D. Poland, S. S. Pufu, D. Simmons-Duffin and R. Yacoby,
  "Bootstrapping 3D Fermions,"
  [arXiv:1508.00012 [hep-th]](http://arxiv.org/abs/1508.00012).

- D. Poland and A. Stergiou,
  "Exploring the Minimal 4D N=1 SCFT,"
  [arXiv:1509.06368 [hep-th]](http://arxiv.org/abs/1509.06368).

-  M. Lemos and P. Liendo,
  "Bootstrapping N=2 chiral correlators,"
  [arXiv:1510.03866 [hep-th]](http://arxiv.org/abs/1510.03866).

(Please let me know if you would like your paper to be included in this list!)
