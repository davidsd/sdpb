# Outer_Limits Usage

The exact details of how `outer_limits` works is described in the
`outer_limits` [manual](Outer_Limits/Outer_Limits.pdf). To briefly summarize,
`outer_limits` reads in function evaluated at Chebyshev zeros and
initial points.  It applies constant constraints at those initial
points and computes a solution.  It then iteratively adds constraints
where the functional is not positive semidefinite and computes new
solutions until the resulting functional is positive semidefinite
everywhere.

## Differences from SDPB

`outer_limits` uses the `sdpb` routines internally when computing the
solution at each iteration, so many of the options are the same as
`sdpb`.  Running `outer_limits --help` will print out a complete list
of options.  There are four main differences:

1. `functions`: Mathematica, JSON, or NSV file with SDP functions evaluated
   at Chebyshev zeros.

2. `points`: JSON or NSV file with initial points.

3. `dualityGapReduction`: Shrink the duality gap threshold by this
   factor during each outer iteration.  Smaller means slower
   convergence.

4. `out`: The optimal solution is saved to this file in json
    format. Defaults to `functions` with the ending `_out.json`.

The format for `functions` is detailed in the
[schema](json_schema/functions_schema.json), and there is an
[example](../test/data/outer_limits/toy_functions.json) included with the source code.
Similarly, `points` has a [schema](json_schema/points_schema.json) and an
[example](../test/data/outer_limits/toy_functions_points.json).

If you have an existing input file for `sdpb`, you can use the bundled
`sdp2functions` and `pvm2functions` programs to create suitable input
from JSON and XML files.  They work very similarly to `sdp2input` and
`pvm2sdp`.  You can run them with the `--help` option to see usage
information.

All of the other options have a similar meaning as for `sdpb`.  For
example, `initialMatrixScaleDual` will set the initial size of the `Y`
matrices when computing a solution for each outer iteration.

## Running Outer_Limits

Running `outer_limits` is very similar to running `sdpb`.  For
example, to run the test case on your quad-core laptop, you can
probably use `mpirun`:

    mpirun -n 4 build/outer_limits --functions test/data/outer_limits/toy_functions.json --points test/data/outer_limits/toy_functions_points.json --precision=64 --dualityGapThreshold=1e-10 --primalErrorThreshold=1e-10 --dualErrorThreshold=1e-10 --initialMatrixScalePrimal=1e1 --initialMatrixScaleDual=1e1 --maxIterations=1000

Note that this requires relatively low precision: `64`.  Also, the
values for `initialMatrixScalePrimal` and `initialMatrixScaleDual` are
close to `1`.  This is in stark contrast to `sdpb`, where the default
for these numbers is `1e20`.  This is enabled by the scaling that
`outer_limits` applies before it solves each inner iteration.
