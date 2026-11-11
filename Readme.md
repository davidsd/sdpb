[![License](https://img.shields.io/github/license/davidsd/sdpb)](LICENSE)
[![Release](https://img.shields.io/github/v/release/davidsd/sdpb)](https://github.com/davidsd/sdpb/releases/latest)
[![CircleCI](https://dl.circleci.com/status-badge/img/gh/davidsd/sdpb/tree/master.svg?style=shield)](https://dl.circleci.com/status-badge/redirect/gh/davidsd/sdpb/tree/master)
[![Issues](https://img.shields.io/github/issues/davidsd/sdpb)](https://github.com/davidsd/sdpb/issues)

## Contents

* [SDPB](#sdpb)
* [Installation and Usage](#installation-and-usage)
* [Experimental branches](#experimental-branches)
* [Test problems](#test-problems)
* [Front ends](#front-ends)
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

SDPB can also solve general semidefinite programs in SDPA format,
see [Solving SDPA problems](docs/Usage.md#solving-sdpa-problems).

Authors: David Simmons-Duffin (dsd@caltech.edu), Walter Landry (wlandry@caltech.edu), Vasiliy Dommes (vasdommes@gmail.com)

## Installation and Usage

The easiest way to run SDPB on a Windows or Mac machine is to follow
the [Docker instructions](docs/Docker.md).  For Linux and HPC centers,
the [Singularity](docs/Singularity.md) instructions will probably work
better.  If you want to build it yourself, there are detailed
instructions in [Install.md](Install.md).

Usage instructions are detailed in [Usage.md](docs/Usage.md).

## Experimental branches

Some new features (faster CPU paths, new options of `sdpb`, `pmp2sdp` and `spectrum`, and a GPU backend) are
offered as pull requests before they are reviewed and merged. [docs/Experimental.md](docs/Experimental.md) lists
them, with their review status and the people to credit if you use them, and says how to build and combine them.

## Test problems

For work on SDPB's performance, there is a ladder of optimization test problems of increasing size,
`sdpb-test-01` to `sdpb-test-08`, from small problems that solve on one node to multi-node and multi-GPU ones.
Each comes with its SDP, its PMP (or a recipe to make it), a reference solution or reference timings, a README
and a check script. They are
[here](https://www.dropbox.com/scl/fo/cvpy6qz6ztsf3ibjpk7k6/AFS2Eb7iwgVNeeUgzm-WNFU?rlkey=0nfp6xddwcbe5creg7q4subiq&st=v7dq9dr5&dl=0);
more rungs are being added, among them SDPA-format problems (`.dat-s` files, which `sdpb` reads directly; see
[Solving SDPA problems](docs/Usage.md#solving-sdpa-problems)). Measured timings on them are collected in
[docs/Performance.md](docs/Performance.md).

## Front ends

Several front ends for SDPB are available:

- [PyCFTBoot](https://github.com/cbehan/pycftboot) (Python) by Connor Behan ([arXiv:1602.02810](https://arxiv.org/abs/1602.02810), [INSPIRE](https://inspirehep.net/literature/1420573))
- [cboot](https://github.com/tohtsky/cboot) (Python) by Tomoki Ohtsuki ([arXiv:1602.07295](https://arxiv.org/abs/1602.07295), [INSPIRE](https://inspirehep.net/literature/1423075))
- [virasoro_blocks](https://github.com/davidpoland/virasoro_blocks) (Python) by David Poland ([arXiv:2610.39106](https://arxiv.org/abs/2610.39106), [INSPIRE](https://inspirehep.net/literature/3220379)), for the 2d bootstrap
- [simpleboot](https://gitlab.com/bootstrapcollaboration/simpleboot) (Mathematica) by Ning Su
- [autoboot](https://github.com/selpoG/autoboot) (Mathematica) by Mocho Go and Yuji Tachikawa ([arXiv:1903.10522](https://arxiv.org/abs/1903.10522), [INSPIRE](https://inspirehep.net/literature/1726773)), which generates bootstrap equations for cboot
- [qboot](https://github.com/selpoG/qboot) (C++) by Mocho Go ([arXiv:2006.04173](https://arxiv.org/abs/2006.04173), [INSPIRE](https://inspirehep.net/literature/1799886))
- [sdpb-haskell](https://gitlab.com/davidsd/sdpb-haskell) (Haskell) by David Simmons-Duffin

## Attribution

### SDPB

If you use SDPB in work that results in publication, consider citing

- D. Simmons-Duffin, *A Semidefinite Program Solver for the Conformal Bootstrap*, JHEP 1506, 174 (2015)
  [arXiv:1502.02033](https://arxiv.org/abs/1502.02033), [INSPIRE](https://inspirehep.net/literature/1343540).
- W. Landry and D. Simmons-Duffin, *Scaling the semidefinite program solver SDPB*
  [arXiv:1909.09745](https://arxiv.org/abs/1909.09745), [INSPIRE](https://inspirehep.net/literature/1755416).

Depending on how SDPB is used, the following might also be relevant:

- First use of semidefinite programming in the bootstrap:
  - D. Poland, D. Simmons-Duffin and A. Vichi,
    *Carving Out the Space of 4D CFTs*, JHEP 1205, 110 (2012) [arXiv:1109.5176](https://arxiv.org/abs/1109.5176), [INSPIRE](https://inspirehep.net/literature/929469).
- Arbitrary spacetime dimension:
  - F. Kos, D. Poland and D. Simmons-Duffin, *Bootstrapping the O(N) Vector Models*,
    JHEP 1406, 091 (2014) [arXiv:1307.6856](https://arxiv.org/abs/1307.6856), [INSPIRE](https://inspirehep.net/literature/1244691).
- Arbitrary systems of correlation functions:
  - F. Kos, D. Poland and D. Simmons-Duffin, *Bootstrapping Mixed
    Correlators in the 3D Ising Model*, JHEP 1411, 109 (2014) [arXiv:1406.4858](https://arxiv.org/abs/1406.4858), [INSPIRE](https://inspirehep.net/literature/1301561).
- Optimal interpolation, for better performance and condition numbers:
  - C.-H. Chang, V. Dommes, P. Kravchuk,
    D. Poland and D. Simmons-Duffin, *Accurate bootstrap bounds from optimal interpolation*, JHEP 2605, 275 (2026)
    [arXiv:2509.14307](https://arxiv.org/abs/2509.14307), [INSPIRE](https://inspirehep.net/literature/2970997).

### approx_objective

- Linear and quadratic variations of the objective function:
  - M. Reehorst, S. Rychkov, D. Simmons-Duffin,
    B. Sirois, N. Su and B. van Rees, *Navigator Function for the Conformal Bootstrap*, SciPost Phys. 11, 072 (2021)
    [arXiv:2104.09518](https://arxiv.org/abs/2104.09518), [INSPIRE](https://inspirehep.net/literature/1859324).

### spectrum

- Original use:
  - Z. Komargodski and D. Simmons-Duffin, *The Random-Bond Ising Model in 2.01 and 3 Dimensions*,
    J. Phys. A 50, 154001 (2017) [arXiv:1603.04444](https://arxiv.org/abs/1603.04444), [INSPIRE](https://inspirehep.net/literature/1427756).
- How it works:
  - D. Simmons-Duffin, *The Lightcone Bootstrap and the Spectrum of the 3d Ising CFT*, JHEP 1703, 086
    (2017) [arXiv:1612.08471](https://arxiv.org/abs/1612.08471), [INSPIRE](https://inspirehep.net/literature/1506449).

## Acknowledgements

- Versions 2 and 3 of SDPB were made possible by support from the [Simons Collaboration on the Nonperturbative Bootstrap](http://bootstrapcollaboration.com/).

- The design of SDPB was partially based on the solvers [SDPA](http://sdpa.sourceforge.net/) and SDPA-GMP, which were essential sources of inspiration and examples.

- Thanks to Filip Kos, David Poland, and Alessandro Vichi for collaboration in developing semidefinite programming methods for the conformal bootstrap and assistance testing SDPB.

- Thanks to Amir Ali Ahmadi, Hande Benson, Pablo Parrilo, and Robert Vanderbei for advice and discussions about semidefinite programming.

- Thanks also to Noah Stein, who first suggested the idea of semidefinite programming to me in [this Math Overflow question](http://mathoverflow.net/questions/33242/continuous-linear-programming-estimating-a-solution).

## Works Using SDPB

Here is [a list of papers citing SDPB](http://inspirehep.net/search?ln=en&p=refersto%3Arecid%3A1343540&sf=earliestdate).
