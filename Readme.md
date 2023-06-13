## Contents

* [Skydive](#skydive)
* [Installation and Usage](#installation-and-usage)
* [Attribution](#attribution)
* [Acknowledgements](#acknowledgements)
* [Works Using SDPB](#works-using-sdpb)

# Skydive

Skydive is an open-source, arbitrary-precision, parallelized solver for the PMP bundle problem: a section of PMP that smoothly depends on a few ``external parameters". The solver attempts to find a point in the space of external parameters such that (1), the objective of PMP is minimized; or (2), the external parameters are extremized along certain direction while keeping the objective of PMP non-positivie.

The description of the strategy and alogrithm is in the paper : 

This code is developed based on SDPB 2.5.1.

## Installation and Usage

The installation of `skydive` is the same as SDPB 2.5.1, except the one need a modified verision of `elemental` library : [https://gitlab.com/AikeLiu/elemental](https://gitlab.com/AikeLiu/elemental/-/tree/MPFR_LogExp) (the MPFR_LogExp branch). After installing the modified `elemental`, follow the rest of instructions in the [sdpb/tree/master/docs/site_installs](https://github.com/davidsd/sdpb/tree/master/docs/site_installs).

This `skydive` program perform one iteration in the main loop of the skydiving algorithm. For conformal bootstrap applications, the main loop is drived and controlled by one of the two external program:
- [Hyperion](https://gitlab.com/davidsd/dynamical-sdp).
- [simpleboot](https://gitlab.com/bootstrapcollaboration/simpleboot).
The external programs set up the bootstrap problem, produce SDPs, call `skydive`, update $p$ and other parameters, and iterate the main loop until the goal is achieved. 

## Attribution

## Acknowledgements



