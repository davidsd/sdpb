# Usage

Details of how SDPB works are described in the [the
manual](SDPB-Manual.pdf).  An example input file
[test.xml](../test/test.xml) is included with the source code.

The build system creates the executables `pvm2sdp`, `sdp2input`, and
`sdpb` in the `build` directory.  There are two steps when running
SDPB.

## Create input files

You will normally start with either an SDP file in Mathematica or JSON format, or
a Polynomial Vector Matrices file in XML format.  These must first be
converted, using `sdp2input` or `pvm2sdp`, into a format that SDPB can
quickly load.  When creating these input files, you must choose a
working precision.  In general, you should use the same precision as
when you run `sdpb`, though you may also use a larger precision.  Both
`sdp2input` and `pvm2sdp` will run faster in parallel.

### Converting an SDP

Use `sdp2input` to create input from files with an SDP.  The usage is

    sdp2input --precision=[PRECISION] --input=[INPUT] --output=[OUTPUT]

`[PRECISION]` is the number of bits of precision used in the
conversion.  `[INPUT]` is a single Mathematica, JSON, or NSV
(Null Separated Value) file.  `[OUTPUT]` is an output directory.

The single file Mathematica and JSON formats are described in Section
3.2 of the [the manual](SDPB-Manual.pdf).  In addition, for JSON there
is a [schema](sdp2input_schema.json).

The NSV format allows you to load an SDP from multiple JSON and/or
Mathematica files.  NSV files contain a list of files, separated by
null's.  When using multiple files, each component is optional.  If
there are multiple versions of `normalization` or `objective`,
`sdp2input` will use the last one it sees.  Multiple elements of
`PositiveMatrixWithPrefactor` will be concatenated together.

One way to generate NSV files is with the `find` command.  For
example, if you have a directory `input` with many files, you can
generate an NSV file with the command

    find input/ -print0 -name "*.m" > file_list.nsv

`sdp2input` assumes that files ending with `.nsv` are NSV,
files ending with `.json` are JSON, and everything else is
Mathematica.  NSV files can also recursively reference other NSV
files.

There are example input files in
[Mathematica](../test/sdp2input_test.m),
[JSON](../test/sdp2input_test.json), and
[NSV](../test/sdp2input_split.nsv) format.  They all define the same
SDP, with the NSV example loading the SDP from two Mathematica files:
[sdp2input\_split1.m](../test/sdp2input_split1.m) and
[sdp2input\_split2.m](../test/sdp2input_split2.m).

### Converting Polynomial Vector Matrices

Use `pvm2sdp` to create input files from Polynomial Vector Matrix files
in XML.  The format for these XML files is described in Section 3.1 of
[the manual](SDPB-Manual.pdf).  The usage is

    pvm2sdp [PRECISION] [INPUT] ... [OUTPUT]

`[PRECISION]` is the number of bits of precision used in the
conversion.  `[INPUT] ...` is a list of one or more files.  Each of
those input files can be either an XML file with Polynomial Vector
Matrices, or an Null Separated Value (NSV) file with a list of
files. `[OUTPUT]` is an output directory.

For example, the command to convert the file `test/test.xml`, using
1024 bits of precision, and store the result in the directory
`test/test/`, is

    pvm2sdp 1024 test/test.xml test/test

For input, you can also use an NSV file containing a list of file names
separated by nulls.  There is an example in `test/file_list.nsv`.
Using it is similar to the xml files

    pvm2sdp 1024 test/file_list.nsv test/test

One way to generate this file list is with the `find` command.  For
example,

    find test -print0 -name "test.xml" > test/file_list.nsv
    
will regenerate `test/file_list.nsv`.  If you have a directory `input`
with many xml files, you can generate a file list with the command

    find input/ -print0 -name "*.xml" > file_list.nsv

`pvm2sdp` assumes that anything with the `.nsv` extension is a null
separated file list, and everything else is an xml file.  File lists
can also recursively reference other file lists.  So running

    find test -print0 -name "*.nsv" > file_list_list.nsv
    pvm2sdp 1024 file_list_list.nsv test/test

will also work.

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
of the time to evaluate each block.  If `block_timings` does not
already exists in the input directory or a checkpoint directory, SDPB
will create one.  SDPB will run for 2 iterations and write the time to
evaluate each block into `block_timings`.  SDPB has to run for 2
iterations because measuring the first step generally gives a poor
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

## Optimizing Memory Use

SDPB's defaults are set for optimal performance.  This may result in
using more memory than is available.  Running SDPB on more nodes will
reduce the amount of memory required on each node.  If this is not
sufficient, you can also also use the option `--procGranularity`.
This option sets minimum number of processes that a block group can
have, so it must evenly divide the `--procsPerNode` option.  Using a
larger granularity will result in less memory use (up to a point)
because SDPB will make fewer local copies of the matrix Q.  However,
larger granularity is also slower because even small blocks will be
distributed among multiple cores.  So you should use
`--procGranularity` only when absolutely needed.
