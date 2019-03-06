# Usage

Details of how SDPB works are described in the [the
manual](SDPB-Manual.pdf).  An example input file
[test.xml](../test/test.xml) is included with the source code.

The build system creates the executables `pvm2sdp`, `sdp2input`, and
`sdpb` in the `build` directory.  There are two steps when running
SDPB.

## Create input files

You will normally start with either a Mathematica file with an SDP, or
an XML file with Polynomial Vector Matrices.  These must first be
converted, using `sdp2input` or `pvm2sdp`, into a format that SDPB can
quickly load.  When creating these input files, you must choose a
working precision.  In general, you should use the same precision as
when you run `sdpb`, though you may also use a larger precision.  Both
`sdp2input` and `pvm2sdp` will run faster in parallel.

### Converting SDP Mathematica files

Use `sdp2input` to create input from Mathematica files with an SDP.
The format for these Mathematica files is described in Section 3.2 of
the [the manual](SDPB-Manual.pdf).  The usage is

    sdp2input --precision=[PRECISION] --input=[INPUT] --output=[OUTPUT]

`[PRECISION]` is the number of bits of precision used in the
conversion.  `[INPUT]` is a Mathematica file with an SDP. `[OUTPUT]`
is an output directory.

### Converting XML files

Use `pvm2sdp` to create input files from XML with polynomial vector
matrices.  The format for these XML files is described in Section 3.1
of [the manual](SDPB-Manual.pdf).  The usage is

    pvm2sdp [PRECISION] [INPUT] ... [OUTPUT]

`[PRECISION]` is the number of bits of precision used in the
conversion.  `[INPUT] ...` is a list of one or more xml
files. `[OUTPUT]` is an output directory.

For example, the command to convert the file `test/test.xml`, using
1024 bits of precision, and store the result in the directory
`test/test/`, is

    pvm2sdp 1024 test/test.xml test/test

## Running SDPB.

The options to SDPB are described in detail in the help text, obtained
by running `build/sdpb --help`.  The most important options are `-s [--sdpDir]`,
`--precision`, and `--procsPerNode`.

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

To efficiently run large MPI jobs, SDPB needs an accurate measurement
of the time to evaluate each block.  If a `block_timings` does not
already exists in the input directory or a checkpoint directory, SDPB
will create one.  SDPB will run for 2 iterations and write the time to
evaluate each block into `block_timings`.  SDPB has to run for 2
iterations because measuring first step generally gives a poor
estimate.  During the first step, many quantities may be zero.
Adding and multiplying zero is much faster with extended precision.

If you are running a large family of input files with the same
structure but different numbers, the measurements are unlikely to
differ.  In that case, you can reuse timings from previous inputs by
copying the `block_timings` file to other input directories.

If different runs have the same block structure, you can also reuse
checkpoints from other inputs. For example, if you have a previous
checkpoint in `test/test.ck`, you can reuse it for a different input
in `test/test2` with a command like

    mpirun -n 4 build/sdpb --precision=1024 --procsPerNode=4 -s test/test2/ -i test/test.ck

In addition to having the same block structure, the runs must also use
the same `precision`, `procsPerNode`, and number and distribution of
cores.
