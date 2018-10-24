# Usage

The input format for SDPB is XML-based and described in
[the manual](SDPB-Manual.pdf).
The Mathematica file [SDPB.m](../mathematica/SDPB.m) includes code to export semidefinite
programs in this format, along with some examples. An example input
file [test.xml](../test/test.xml) is included with the source code.

The build system creates the executables `pvm2sdp` and `sdpb` in the
build directory.  There are three steps when running SDPB.

## Converting XML files

XML files with polynomial vector matrices must first be converted with
`pvm2sdp` into a format that SDPB can quickly load.  The usage is

    build/pvm2sdp [PRECISION] [INPUT] ... [OUTPUT]

`[PRECISION]` is the number of bits of precision used in the
conversion.  `[INPUT] ...` is a list of one or more xml
files. `[OUTPUT]` is an output directory.

For example, the command to convert the file `test/test.xml`, using
1024 bits of precision, and store the result in the directory
`test/test/`, is

    build/pvm2sdp 1024 test/test.xml test/test
  
`pvm2sdp` writes the output in a text format, so it is perfectly OK to
use a larger precision with `pvm2sdp` than when running `sdpb`.

## [Optional] Measure time to evaluate each block.

To efficiently run large MPI jobs, SDPB needs an accurate measurement
of the time to evaluate each block.  You may skip this step if you
have many, many more blocks than cores.

To measure the timings, run SDPB on your input for two steps.  Two
steps is required because the first step may have very different
performance characteristics.  During the first step, many quantities
may be zero.  Adding and multiplying zero is much faster with extended
precision.

For the measurements, you should run with the number of cores less
than or equal to the number of blocks.  Otherwise, SDPB will attempt
to split the blocks between cores, which will result in inaccurate
timings.  Once you have timings, you can then run SDPB on as many
cores as you like.

This also means that the measurement run should generally not
checkpoint the results.  Checkpoints are not usable when running
with a different number of cores.

For accurate timings, you should use the same precision as you will
use for the final runs.

As an example, the command to measure timings for `test/test` using a
single core and store the result in `test/test/block_timings` is

    build/sdpb --precision=1024 --maxIterations=2 --noFinalCheckpoint --procsPerNode=1 -s test/test/ --writeBlockTiming

The test input only has a single block, so there should be a single number
in `test/test/block_timings`.

If you are running a large family of input files with the same
structure but different numbers, the measurements are unlikely to
differ.  In that case, you can reuse timings from previous inputs by
copying the `block_timings` file to other input directories.

## Running SDPB.

The options to SDPB are described in detail in the help text, obtained
by running `build/sdpb --help`.  The most important options are `-s [--sdpDir]`,
`--precision`, and `--procsPerNode`.

When using checkpoints, you can only load checkpoints that use the
same `precision`, `procsPerNode`, and number of cores.

SDPB uses MPI to run in parallel, so you may need a special syntax to
launch it.  For example, if you compiled the code on your own laptop,
you will probably use `mpirun` to invoke SDPB.  If you have 4 physical
cores on your machine, the command is

    mpirun -n 4 build/sdpb --precision=1024 --procsPerNode=4 -s test/test/

On the Yale Grace cluster, the command used in the Slurm batch file is

    mpirun build/sdpb --precision=1024 --procsPerNode=$SLURM_NTASKS_PER_NODE -s test/test/

In contrast, the Harvard Odyssey 3 cluster, which also uses Slurm,
uses the srun command

    srun -n $SLURM_NTASKS --mpi=pmi2 build/sdpb --precision=1024 --procsPerNode=$SLURM_NTASKS_PER_NODE -s test/test

The documentation for your HPC system will tell you how to write a
batch script and invoke MPI programs.

