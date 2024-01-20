## Contents

* [Skydive](#skydive)
* [Installation and Usage](#installation-and-usage)

# Skydive

Skydive is an open-source, arbitrary-precision, parallelized solver for the PMP bundle problem: a section of PMP that
smoothly depends on a few "external parameters". The solver attempts to find a point in the space of external parameters
such that (1), the objective of PMP is minimized; or (2), the external parameters are extremized along certain direction
while keeping the objective of PMP non-positive.

The description of the strategy and algorithm is in the paper :

This code is developed based on SDPB 2.5.1 and then updated to SDPB 2.6.1.

## Installation and Usage

The installation of `skydive` is the same as for SDPB 2.6.1,
see instructions at [docs/site_installs](https://github.com/davidsd/sdpb/tree/master/docs/site_installs).

NB: make sure that you are using an up-to-date version
of [Elemental](https://gitlab.com/bootstrapcollaboration/elemental) library
that [implements](https://gitlab.com/bootstrapcollaboration/elemental/-/merge_requests/4)  `Exp(BigFloat)`
and `Log(BigFloat)` functions.

This `skydive` program perform one iteration in the main loop of the skydiving algorithm. For conformal bootstrap
applications, the main loop is driven and controlled by one of the two external program:
- [Hyperion](https://gitlab.com/davidsd/dynamical-sdp).
- [simpleboot](https://gitlab.com/bootstrapcollaboration/simpleboot).

The external programs set up the bootstrap problem, produce SDPs, call `skydive`, update `p` and other parameters, and
iterate the main loop until the goal is achieved. 
