[![License](https://img.shields.io/github/license/davidsd/sdpb)](LICENSE)
[![Release](https://img.shields.io/github/v/release/davidsd/sdpb)](https://github.com/davidsd/sdpb/releases/latest)
[![CircleCI](https://dl.circleci.com/status-badge/img/gh/davidsd/sdpb/tree/master.svg?style=shield)](https://dl.circleci.com/status-badge/redirect/gh/davidsd/sdpb/tree/master)
[![Issues](https://img.shields.io/github/issues/davidsd/sdpb)](https://github.com/davidsd/sdpb/issues)

## Contents

* [SDPB](#sdpb)
* [Installation and Usage](#installation-and-usage)
* [Attribution](#attribution)
* [Acknowledgements](#acknowledgements)
* [Works Using SDPB](#works-using-sdpb)

# SDPB

SDPB is an open-source, arbitrary-precision, parallelized semidefinite
program solver, designed for the conformal bootstrap. It solves the following problem:

Let $S^{m\times m}[x]$ be the space of symmetric $m\times m$ matrices whose entries are polynomials in $x$.

- Given:
  - A collection of $J\times (N+1)$ polynomial matrices $M_{nj}(x) \in S^{m_j\times m_j}[x]$, ($j=1,\dots,J$ and $n=0,\dots,N$),
  - A constant $b_0 \in \mathbb{R}$,
  - A vector $b\in \mathbb{R}^N$,
- Maximize $b_0 + b\cdot y$ over decision variables $y\in \mathbb{R}^N$
- Such that: $M_{0j}(x) + \Sigma_{n=1}^N y_n M_{nj}(x) \succeq 0$ for all $x\geq 0$ and $j=1,\dots,J$.

Here, $M\succeq 0$ means "M is positive semidefinite."

For more information, see [A Semidefinite Program Solver for the Conformal Bootstrap](http://arxiv.org/abs/1502.02033)
and [the manual](/docs/SDPB_Manual/SDPB-Manual.pdf).

Authors: David Simmons-Duffin (dsd@caltech.edu), Walter Landry (wlandry@caltech.edu).

On April 25, 2019, the main branch of this repository was updated to SDPB 2 which has very different performance characteristics and somewhat different usage instructions from version 1.0. Please see the changelog and other documentation for details.

## Installation and Usage

The easiest way to run SDPB on a Windows or Mac machine is to follow
the [Docker instructions](docs/Docker.md).  For Linux and HPC centers,
the [Singularity](docs/Singularity.md) instructions will probably work
better.  If you want to build it yourself, there are detailed
instructions in [Install.md](Install.md).

Usage instructions are detailed in [Usage.md](docs/Usage.md).

Two python wrappers for SDPB are available:

- [PyCFTBoot](https://github.com/cbehan/pycftboot) by Connor Behan ([arXiv:1602.02810](http://arxiv.org/abs/arXiv:1602.02810))
- [cboot](https://github.com/tohtsky/cboot) by Tomoki Ohtsuki ([arXiv:1602.07295](http://arxiv.org/abs/arXiv:1602.07295)).

An unofficial Haskell wrapper is available:

- [sdpb-haskell](https://gitlab.com/davidsd/sdpb-haskell) by David Simmons-Duffin

## Attribution

### SDPB

If you use SDPB in work that results in publication, consider citing

- D. Simmons-Duffin, *A Semidefinite Program Solver for the
  Conformal Bootstrap*, JHEP 1506, 174 (2015) [arXiv:1502.02033](http://arxiv.org/abs/1502.02033).
- W. Landry and D. Simmons-Duffin, *Scaling the semidefinite program solver SDPB*
  [arXiv:1909.09745](https://arxiv.org/abs/1909.09745).

Depending on how SDPB is used, the following other sources might be relevant:

The first use of semidefinite programming in the bootstrap:

- D. Poland, D. Simmons-Duffin and A. Vichi, *Carving Out the Space of
  4D CFTs*, JHEP 1205, 110 (2012) [arXiv:1109.5176](http://arxiv.org/abs/1109.5176).

The generalization of semidefinite programming methods to arbitrary
spacetime dimension:

- F. Kos, D. Poland and D. Simmons-Duffin, *Bootstrapping the O(N)
  Vector Models*, JHEP 1406, 091 (2014) [arXiv:1307.6856](http://arxiv.org/abs/1307.6856).

The generalization of semidefinite programming methods to arbitrary
systems of correlation functions:

- F. Kos, D. Poland and D. Simmons-Duffin, *Bootstrapping Mixed
  Correlators in the 3D Ising Model*, JHEP 1411, 109 (2014) [arXiv:1406.4858](http://arxiv.org/abs/1406.4858).

### approx_objective

Derivation of linear and quadratic variations of the objective function, used in `approx_objective`:

- M. Reehorst, S. Rychkov, D. Simmons-Duffin, B. Sirois, N. Su, B. van Rees, *Navigator Function for the Conformal Bootstrap*,
  [arXiv:2104.09518](http://arxiv.org/abs/2104.09518).

### spectrum

Spectrum extraction was originally written for use in:

  - Z. Komargodski and D. Simmons-Duffin, *The Random Bond
    Ising Model in 2.01 and 3 Dimensions*, [arXiv:1603.04444](https://arxiv.org/abs/1603.04444)

An explanation of how it works appears in:

  - D. Simmons-Duffin, *The Lightcone Bootstrap and the Spectrum of the 3d Ising CFT*, [arXiv:1612.08471](https://arxiv.org/abs/1612.08471)

## Acknowledgements

- Version 2 of SDPB was made possible by the [Simons Collaboration on the Nonperturbative Bootstrap](http://bootstrapcollaboration.com/).

- The design of SDPB was partially based on the solvers [SDPA](http://sdpa.sourceforge.net/) and SDPA-GMP, which were essential sources of inspiration and examples.

- Thanks to Filip Kos, David Poland, and Alessandro Vichi for collaboration in developing semidefinite programming methods for the conformal bootstrap and assistance testing SDPB.

- Thanks to Amir Ali Ahmadi, Hande Benson, Pablo Parrilo, and Robert Vanderbei for advice and discussions about semidefinite programming.

- Thanks also to Noah Stein, who first suggested the idea of semidefinite programming to me in [this Math Overflow question](http://mathoverflow.net/questions/33242/continuous-linear-programming-estimating-a-solution).

## Works Using SDPB

As of April 2019, SDPB has been used in approximately 70 works. Here is [a list of papers citing SDPB](http://inspirehep.net/search?ln=en&p=refersto%3Arecid%3A1343540&sf=earliestdate).
