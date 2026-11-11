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

`sdpb` can also solve general semidefinite programs in SDPA format, see [Solving SDPA problems](#solving-sdpa-problems).

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
You can also use `--maxNumPoles` to downsample the PMP (see Section 4.3 in https://arxiv.org/abs/2509.14307), and `--outputFormat` and `--zip` to specify output file format. Call `pmp2sdp --help` to see options description. 

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

## Solving SDPA problems

`sdpb` can also solve general (non-polynomial) semidefinite programs in the SDPA format,
which is used e.g. by [SDPA](http://sdpa.sourceforge.net/), [SDPA-GMP](https://github.com/nakatamaho/sdpa-gmp)
and the [SDPLIB](../test/data/SDPLIB/README.md) collection of test problems.
Pass the file to `-s [--sdpDir]` instead of the `sdp` directory created by `pmp2sdp`:

    mpirun -n 4 build/sdpb --precision=512 -s problem.dat-s -o out -c ck

The file extension selects the format:

- `.dat-s`: sparse format, see [SDPA sparse format](../test/data/SDPLIB/README.md#sdpa-sparse-format).
- `.dat`: dense format, where each matrix is written in full, e.g. `{ {1, 2}, {2, 3} }`.
- `.dat-s.gz`, `.dat.gz`: the same formats, compressed with gzip.

As in SDPA, the problem is:

- Primal: minimize $c_0 + \sum_{i=1}^m c_i x_i$ over $x \in \mathbb{R}^m$, such that $X = \sum_{i=1}^m F_i x_i - F_0 \succeq 0$.
- Dual: maximize $c_0 + F_0 \bullet Y$ over $Y \succeq 0$, such that $F_i \bullet Y = c_i$ for $i=1,\dots,m$.

The matrices $F_i$ are block diagonal. A negative block size in the input denotes a diagonal block;
SDPB splits it into $1\times 1$ blocks.

**Constant objective term.** The standard SDPA format has no constant term, so $c_0=0$ by default.
SDPB adds an extension: to set $c_0$, add a comment line `"sdpb.objectiveConst <value>"`
(or `* sdpb.objectiveConst <value>`) among the comment lines at the beginning of the file, e.g.

    "Example 1: mDim = 3, nBLOCK = 1, {2}"
    "sdpb.objectiveConst -0.1"
       3  =  mDIM
       ...

Other SDPA solvers ignore this line, since it is a comment.
See [example1.dat](../test/data/end-to-end_tests/sdpa/example1-objective-const/input/example1.dat).

**Output.** `out.txt` has the same fields as for SDPB problems.
`primalObjective` is $c_0 + \sum_i c_i x_i$, and `dualObjective` is $c_0 + F_0 \bullet Y$.
The `--writeSolution` option accepts:

- `x`: the vector $x$, written to `x.txt`.
- `X`, `Y`: the matrices $X$ and $Y$, one file per block: `X_matrix_0.txt`, `X_matrix_1.txt`, etc.
  Block indices start from 0. A diagonal block is split into $1\times 1$ blocks, and each of them is written to a separate file.

The default is `--writeSolution=x`. `y` and `z` are not available for SDPA problems.

**Other options.** All solver options (`--precision`, thresholds, `--initialMatrixScalePrimal`, `--maxIterations`,
`--maxMemory` etc.) have the same meaning as for SDPB problems. Checkpoints work in the same way.
SDPB reads and writes `block_timings` only in the checkpoint directory, since the input is a file.
`approx_objective` and `spectrum` work only with SDPB problems created by `pmp2sdp`.

**Performance.** The SDPA solver is new and currently slower than
[sdpa-gmp-omp](https://github.com/Canonical111/sdpa-gmp-omp) (a fork of SDPA-GMP with fixed OpenMP parallelization):
in [their benchmark](https://github.com/Canonical111/sdpa-gmp-omp/blob/main/BENCHMARKS.md#versus-sdpb-the-standard-multiprecision-sdp-solver),
SDPB needed 2.4x more iterations and about 2x more time per iteration, and used 25x more RAM.
Until this is fixed, we recommend sdpa-gmp-omp for all problems that fit on a single node.

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

    mpirun -n 4 build/spectrum --pmpInfo=test/data/end-to-end_tests/1d/output/sdp/pmp_info.json --output=test/out/spectrum/1d/pmp.json/spectrum.json --precision=768 --solution=test/data/end-to-end_tests/1d/output/out --threshold=1e-10

This will output the spectra into `test/out/spectrum/1d/spectrum.json` and should look like

```
[
  {
    "block_path": "test/data/end-to-end_tests/1d/input/pmp.json",
    "zeros": [
      {
        "zero": "1.04249678571815812098400651940402588551817464537296262014475782848711468517831729605366065940294278580341286805477340769648758871722958945799926776073918390590815495720242896470413021851113098588836134728703346283",
        "lambda": [
          [
            "1.54694357833877357195864820903901116949337414099884811266823030281476233526526217385572091950844480308078719637826381700862165648847084444755689823793327877973940615930065072306902959073017284495666199126761059837"
          ]
        ]
      }
    ],
    "error": "2.61558507636706998768268589785574818292702184050960721508531240370407569660416817267708593353995302214421134564309862375671483311967442639952542078221319922448465353968305149652407922145661156680520217679074769755e-26"
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

`"lambda"` in `spectrum.json` is an array of OPE coefficient vectors
(eigenvectors of $V_{j,\tau}$ with non-negligible eigenvalues, see footnote 47 in [arxiv:1612.08471](https://arxiv.org/abs/1612.08471)),
defined as
```math
\vec{\lambda}_{j,x} = \vec{v}_{j,x} / \sqrt{\chi_j^\prime(x)}
```
where $\vec{v}_{j,x}$ is given by Eq. (A.7) in [arxiv:1612.08471](https://arxiv.org/abs/1612.08471)
and $\chi_j^\prime(x)$ is a `reducedPrefactor` from pmp.json (see [SDPB Manual](SDPB_Manual/SDPB-Manual.pdf) for PMP format description).
Note that this definition disagrees with Eq. (A.8) in [arxiv:1612.08471](https://arxiv.org/abs/1612.08471), which is incorrect.

Due to numerical errors, $V_{j,\tau}$ can have small non-zero eigenvalues `~ O(dualityGap) * max_eigenvalue`.
To filter them out, you may use the option `--minEigenvalueRatio`.
By default, `spectrum` will set `--minEigenvalueRatio` to `sqrt(dualityGap)`.

`spectrum` detects zeros by checking conditions like `f(x) / f(x_1) < threshold`.
You may adjust its sensitivity by changing the option `--threshold`. By default, it is set to `sqrt(dualityGap)`.

If you get several zeros close to each other and want to replace them with a single zero, increase `--minZeroDistance`.

If `spectrum` finds a large unphysical zero, set `--maxZero` to a lower value to exclude it.

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

How to reduce memory usage:

1. Running SDPB on more nodes will reduce the amount of memory required on each node.
2. Decrease `--maxMemory`, e.g. `--maxMemory=60%` or `--maxMemory=64G`. This is the limit for the estimated memory usage on each node. By default, it is 80% of `MemAvailable` at program start. To stay within the limit, SDPB splits shared memory windows used for matrix multiplication, see [bigint_syrk/Readme.md](../src/sdp_solve/SDP_Solver/run/bigint_syrk/Readme.md) for details. If even maximal splitting is not enough, SDPB stops with a "Not enough memory" error.
3. You can also limit only the shared memory windows with `--maxSharedMemory`, e.g. `--maxSharedMemory=32G`. By default, there is no separate limit for them.
4. If some large SDP blocks are assigned to too many cores, this can lead to a significant memory overhead. The worst case is when the number of cores is a huge prime number (e.g. 13 or 17). In that scenario, matrices are distributed across a 1D MPI grid (e.g. `13x1` or `17x1`), which is generally bad for performance and memory (in particular, for the function `El::Trsm()`). Set e.g. `--procGranularity 2` (or `--procGranularity 4`) and/or reduce the number of cores per node. To see block distribution, run SDPB with `--verbosity debug`.

Note that memory estimates are imprecise, and actual memory usage can be higher than expected. If the default `--maxMemory` does not prevent OOM, consider decreasing it and/or increasing number of nodes.
To print detailed memory estimates, run SDPB with `--verbosity debug`.

Smaller limits may affect performance. If shared memory windows have to be split, SDPB prints a warning with current and optimal window sizes.
In our benchmarks, the negative effect on performance was significant only when the shared memory limit was much smaller than the output window size.

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
