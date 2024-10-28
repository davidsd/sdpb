# Usage

Details of how SDPB works are described in the
[manual](SDPB_Manual/SDPB-Manual.pdf). An example input file
[pmp.json](../test/data/end-to-end_tests/1d/input/pmp.json) is included with the source code.

Some known issues and workaround are described [below](#common-issues-and-workarounds).
You may also [find](https://github.com/davidsd/sdpb/issues) unresolved issues
or [report](https://github.com/davidsd/sdpb/issues/new) a new one in the GitHub repository.

The build system creates the executables `pmp2sdp` and
`sdpb` in the `build` directory.  There are two steps when running
SDPB.

## Create input files

You will normally start with a Polynomial Matrix Program (PMP) described in a file (or several files) in JSON,
Mathematica, or XML format. These
files
must first be converted, using `pmp2sdp`, into an SDP format
that SDPB can quickly load.  The format is described in
[SDPB_input_format.md](SDPB_input_format.md).  When creating these
input files, you must choose a working precision. You
should use the same precision as when you run `sdpb` (you may also use larger precision if you are using `pmp2sdp`
with `--outputFormat=json`). `pmp2sdp` will run
faster in parallel.

### Converting PMP to SDP

Use `pmp2sdp` to create SDPB input from files with a PMP. The usage is

    pmp2sdp --precision=[PRECISION] --input=[INPUT] --output=[OUTPUT]

`[PRECISION]` is the number of bits of precision used in the
conversion.  `[INPUT]` is a single Mathematica, JSON, XML or NSV
(Null Separated Value) file. `[OUTPUT]` is an output directory.

The single file Mathematica and JSON formats are described in Section
3.2 of the [manual](SDPB_Manual/SDPB-Manual.pdf). In addition, for JSON there
is a [schema](json_schema/pmp_schema.json).
The format for the XML files is described in Section 3.1 of
[the manual](SDPB_Manual/SDPB-Manual.pdf).

The NSV format allows you to load a PMP from multiple files.
NSV files contain a list of files, separated by null's (`'\0'`).
When using multiple files, each component is optional.
Multiple elements of `PositiveMatrixWithPrefactor` will be concatenated together.
Multiple versions of `normalization` or `objective` are allowed as long as they are identical.

One way to generate NSV files is with the `find` command.  For
example, if you have a directory `input` with many files, you can
generate an NSV file with the command

    find input/ -name "*.m" -print0 > file_list.nsv

`pmp2sdp` assumes that files ending with `.nsv` are NSV,
files ending with `.json` are JSON, `.m` is
Mathematica, and `.xml` is XML.
NSV files can also recursively reference other NSV files.

There is an example [pmp.json](../test/data/end-to-end_tests/1d/input/pmp.json) with a simple one-dimensional PMP
described in manual.
Other JSON files in the folder [test/data/end-to-end_tests/1d/input/](../test/data/end-to-end_tests/1d/input/)
illustrate different ways to define the same PMP.

There are also example input files in
[Mathematica](../test/data/pmp2sdp/m/pmp.m),
[JSON](../test/data/pmp2sdp/json/pmp.json), and
[NSV](../test/data/pmp2sdp/m/file_list.nsv) format. They all define the same
SDP (having three blocks), with the NSV example loading the PMP from two Mathematica files:
[pmp\_split1.m](../test/data/pmp2sdp/m/pmp_split1.m) and
[pmp\_split2.m](../test/data/pmp2sdp/m/pmp_split2.m).

## Running SDPB.

The options to SDPB are described in detail in the help text, obtained
by running `build/sdpb --help`. The most important options are `-s [--sdpDir]` and
`--precision`.
You can specify output and checkpoint directories by `-o [ --outDir ]` and `-c [ --checkpointDir ]`, respectively.

SDPB uses MPI to run in parallel, so you may need a special syntax to
launch it.  For example, if you compiled the code on your own laptop,
you will probably use `mpirun` to invoke SDPB.  If you have 4 physical
cores on your machine, the command is

    mpirun -n 4 build/sdpb --precision=1024 -s test/data/sdp.zip -o test/out/sdpb -c test/out/sdpb/ck

On the Yale Grace cluster, the command used in the Slurm batch file is

    mpirun build/sdpb --precision=1024 -s test/data/sdp.zip -o test/out/sdpb -c test/out/sdpb/ck

In contrast, the Harvard Odyssey 3 cluster, which also uses Slurm,
uses the srun command

    srun -n $SLURM_NTASKS --mpi=pmi2 build/sdpb --precision=1024 -s test/data/sdp.zip -o test/out/sdpb -c test/out/sdpb/ck

The documentation for your HPC system will tell you how to write a
batch script and invoke MPI programs.
Note also that usually you have to **load modules** on your HPC before running SDPB. You can find the corresponding
command in installations instructions for your HPC (see [docs/site_installs](site_installs) folder). For example,
on [Expanse](site_installs/Expanse.md) this command reads

    module load cpu/0.15.4 gcc/10.2.0 openmpi/4.0.4 gmp/6.1.2 mpfr/4.0.2 cmake/3.18.2 openblas/dynamic/0.3.7


Note that most computation for different blocks can be done in parallel, and optimal performance is generally achieved
when the number of MPI jobs is comparable to the number of blocks.

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
checkpoint in `test/out/test.ck`, you can reuse it for a different input
in `test/data/sdp2.zip` with a command like

    mpirun -n 4 build/sdpb --precision=1024 -s test/data/sdp2.zip -i test/out/test.ck

In addition to having the same block structure, the runs must also use
the same `precision`, and number and distribution of cores.

## Running approx_objective

If you have a family of SDP's and a solution to one of these SDP's,
`approx_objective` can compute the approximate value of the objective
for the rest of the family.  The approximation is, by default,
quadratic in the difference between the two SDP's (`b`, `c`, `B`).
`approx_objective` assumes that the bilinear bases `A` are the same.

To compute approximate objectives, write out a text checkpoint when
computing the initial solution with `sdpb` by including the option
`--writeSolution=x,y,X,Y`.  `approx_objective` will then read in this
solution, setup a solver, and compute the new objective.

You specify the location of the new SDP with the option `--newSdp`.
This file would be the output of `pmp2sdp`.
If you have multiple SDP's, then you should create an NSV as in the instructions
for `pmp2sdp` that list all the files.

Setting up the solver can take a long time.  If you have an SDP that
you have perturbed in many different directions, and for logistical
reasons you need to run `approx_objective` separately for each one,
you can save time by saving the solver state with the option
`--writeSolverState`.  To that end, you can also run
`approx_objective` without any new SDP's, only saving the solver
state.

A full example of the whole sequence is

    mpirun -n 4 build/sdpb --precision=1024 -s test/data/sdp -o test/out/approx_objective --writeSolution=x,y,X,Y
    mpirun -n 4 build/approx_objective --precision=1024 --sdp test/data/sdp --writeSolverState
    mpirun -n 4 build/approx_objective --precision=1024 --sdp test/data/sdp --newSdp=test/data/sdp2

The output is a JSON list with each element including the location of
the new SDP, the approximate objective, and the first and second order
terms. There is a [JSON schema](json_schema/approx_objective_schema.json)
describing the format.

The objective can be surprisingly sensitive to small changes,
so the first and second order terms should give you an idea of how
accurate the approximation is.  In general, the first order term
`d_objective`, should be smaller than `objective`.  Also, the second
order term `dd_objective` should be similar in size to the square
of `d_objective` scaled by `objective`.

    dd_objective ~ (d_objective / objective)²

If this is not the case (e.g. `dd_objective > d_objective`), then the
approximate objective is probably inaccurate.

If the linear approximation is accurate enough for you, you can use
the `--linear` option.  This avoids the one-time cost of an expensive
solve.  Also, it only requires the solutions for `x` and `y`, which
SDPB writes out by default.

## Computing the spectrum

Another thing that you can do now that you have a solution is to extract the spectrum.
As a simple example, extracting the spectrum from the toy example would be

    mpirun -n 4 build/spectrum --input=test/data/end-to-end_tests/1d/input/pmp.json --output=test/out/spectrum/1d/spectrum.json --precision=768 --solution=test/data/end-to-end_tests/1d/pmp.json/out --threshold=1e-10

This will output the spectra into `test/out/spectrum/1d/spectrum.json` and should look like

```
[
  {
    "block_path": "test/data/end-to-end_tests/1d/input/pmp.json",
    "zeros":
      [
        {
          "zero": "1.0424967857181581209840065194040256159993020360900482727878770557146614245818844773397565119010581109698382853498679936546513923307745546360686597437951864152447392489871552675002522402284639707590108004747402926563347806805408627031",
          "lambda":
            [
              "1.54694357833877357195864820903901085838088924863264895636292566812483741243475458164155073771873903702617821828504900956715888510611718238011758262357999879234060791367950657551753978299255767817752180863347673614"
            ]
        }
      ],
    "error": "2.6155851084748106058372014479417985936837231505075543152693720913161268299244324614060365156178452537195052377426866117722568851267463995166401153416001203798485032869542329200019745454151278916593515227121025871432270728599938064247e-26"
  }
]
```
    

It is a json file with arrays of zeros. There is a [JSON schema](json_schema/spectrum_schema.json)
describing the format.

The options are described in more detail in the
help text, obtained by running `spectrum --help`.

The spectrum extraction algorithm is described in
[arxiv:1612.08471](https://arxiv.org/abs/1612.08471) (see Appendix A)
and originally implemented in Python, see https://gitlab.com/bootstrapcollaboration/spectrum-extraction.

The vector `"lambda"` in `spectrum.json` is defined as
```math
\vec{\lambda}_{j,x} = \vec{v}_{j,x} / \sqrt{\chi_j^\prime(x)}
```
where $\vec{v}_{j,x}$ is given by Eq. (A.7) in [arxiv:1612.08471](https://arxiv.org/abs/1612.08471)
and $\chi_j^\prime(x)$ is a `reducedPrefactor` from pmp.json (see [SDPB Manual](SDPB_Manual/SDPB-Manual.pdf) for PMP format description).
Note that this definition disagrees with Eq. (A.8) in [arxiv:1612.08471](https://arxiv.org/abs/1612.08471), which is incorrect.

## Common issues and workarounds

### SDPB is slow, how many cores should I use for optimal performance?

Most computation for different blocks can be done in parallel, and optimal performance is generally achieved when the
number of MPI jobs approaches the number of blocks.

Note, however, that increasing number of MPI processes increases also communication overhead, especially between
different machines. Thus, sometimes single-node computation can outperform multi-node ones.

You may use these considerations as a starting point, and run benchmarks in your environment to find the best
configuration for your problem.

### SDPB fails with out-of-memory, std::bad_alloc etc.

SDPB's defaults are set for optimal performance. This may result in using more memory than is available.

Two ways to reduce memory usage:

1. Running SDPB on more nodes will reduce the amount of memory required on each node.
2. Set `--maxSharedMemory` option, e.g. `--maxSharedMemory=64G`. This will reduce memory usage by splitting shared memory windows used for matrix multiplication, see [bigint_syrk/Readme.md](../src/sdp_solve/SDP_Solver/run/bigint_syrk/Readme.md) for details.

If `--maxSharedMemory` is not set by user, SDPB will calculate it automatically based on expected memory usage and amount of available RAM (search for `--maxSharedMemory` in the output to see the new limit).
Note that these estimates are very imprecise, and actual memory usage can be much higher than expected. If automatically calculated `--maxSharedMemory` value does not prevent OOM, consider decreasing it manually and/or increasing number of nodes.

Decreasing `--maxSharedMemory` may affect performance. If the value is too small, SDPB will print a warning with current and optimal shared windows sizes.
In our benchmarks, the negative effect on performance was significant only when `--maxSharedMemory` was much smaller than the output window size.

### SDPB crashes when using all available cores on the node

For older SDPB versions (e.g. 2.5.1), we sometimes observed unexpected crashes for large SDPB runs even with enough memory, e.g. using all 128 cores per node on Expanse
HPC.
In such cases, reducing `$SLURM_NTASKS_PER_NODE` (if you are using SLURM) e.g. from 128 to 64 may help.

### SDPB fails to read large sdp.zip

Sometimes this happens if sdp.zip size exceeds 4GB. You may either regenerate sdp with `pmp2sdp` without `--zip` flag, or unzip sdp manually and pass the resulting folder to `sdpb`:

```
unzip -o path/to/sdp.zip -d path/to/sdp_dir
sdpb -s path/to/sdp_dir <...>
```

### 'A was not numerically HPD' error

Elemental throws this error if you are trying to compute Cholesky decomposition of a matrix that is not Hermitian positive definite (HPD).
Try increasing SDPB `--precision` and/or precision of your PMP input files.

### 'Illegal instruction' crash

This crash has been observed when FLINT binaries compiled on one CPU are used on another CPU which does not support some extensions (e.g. AVX).
For example, this may happen in HPC environment when the code is compiled on a login node and used on a compute node.

Workaround: rebuild FLINT, prodiving `--host` option for `configure` script, e.g. `./configure --host=amd64 <...>`.

See more details in https://github.com/davidsd/sdpb/issues/235.

### 'Running out of inodes' error

This may happen on some filesystems if SDP contains too many block files. Run `pmp2sdp` with `--zip` flag to put
everything into a single zip archive.

### OpenBLAS fails to create threads

```
OpenBLAS blas_thread_init: pthread_create failed for thread 15 of 32: Resource temporarily unavailable
```

Each BLAS call in SDPB should be single-threaded.
To ensure that, disable threading by setting environment variable
```
export OPENBLAS_NUM_THREADS=1
```
and/or
```
export OMP_NUM_THREADS=1
```

### Spectrum does not find zeros

Try to set `--threshold` option for `spectrum` larger than `--dualityGapThreshold` for `sdpb`.

Note that currently spectrum [cannot find isolated zeros](https://github.com/davidsd/sdpb/issues/153).
