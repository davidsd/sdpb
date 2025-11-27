# System info
    Red Hat Enterprise Linux release 8.8 (Ootpa)
    Linux login1.grace.ycrc.yale.edu 4.18.0-477.36.1.el8_8.x86_64 #1 SMP Thu Nov 9 08:12:18 EST 2023 x86_64 x86_64 x86_64 GNU/Linux

### Useful links:

- [Yale HPC documentation](https://docs.ycrc.yale.edu/clusters-at-yale/)
- [Run Jobs with Slurm](https://docs.ycrc.yale.edu/clusters-at-yale/job-scheduling/)
- [Applications & Software ](https://docs.ycrc.yale.edu/clusters-at-yale/applications/)

# Use existing SDPB installation

## Load modules

    module use /home/vd289/install/modules
    module load sdpb

You may also add this command to your `~/.bashrc` file, so that modules will load automatically.

You may run `module -t list` to view loaded modules, and `module purge` to unload all modules.

To see information about the currently loaded module, run

    module whatis sdpb

### Choose SDPB version

By default, `module load sdpb` is equivalent to `module load sdpb/master`.
It loads SDPB built from the latest [master](https://github.com/davidsd/sdpb/tree/master) branch (
run `sdpb --version` to see commit hash, e.g. `SDPB 3.0.0-171-gc39fd506`),
You may list all available versions via

    module av sdpb

and load a specific version via `module load sdpb/<VERSION>`, for example

    module load sdpb/master    
    module load sdpb/3.0.0

## Run SDPB

An example script can be found at `$SDPB_HOME/share/sdpb_example.sh`.
The command

    sbatch $SDPB_HOME/share/sdpb_example.sh

submits the script to the [queueing system](https://docs.ycrc.yale.edu/clusters-at-yale/job-scheduling/).

`sdpb_example.sh` loads modules and runs `pmp2sdp`+`sdpb` for a simple problem.
See script code and comments for more details.

Script output is written to the log file in the current directory, e.g.:
`./sdpb_example.sh.26306151.out`.
SDPB output files are written to the `./out/` folder in the current directory.

# Build SDPB from sources

## Load modules

For compiling and/or running your custom SDPB installation, you have to load modules first:

    module load Autoconf/2.71-GCCcore-12.2.0 Automake/1.16.5-GCCcore-12.2.0 Bison/3.8.2-GCCcore-12.2.0 Boost/1.83.0-GCC-12.2.0 CMake/3.24.3-GCCcore-12.2.0 GMP/6.2.1-GCCcore-12.2.0 libarchive/3.6.1-GCCcore-12.2.0 libtool/2.4.7-GCCcore-12.2.0 MPFR/4.2.0-GCCcore-12.2.0 OpenBLAS/0.3.21-GCC-12.2.0 OpenMPI/4.1.4-GCC-12.2.0 pkg-config/0.29.2-GCCcore-12.2.0 Python/3.10.8-GCCcore-12.2.0 RapidJSON/1.1.0-GCCcore-12.2.0

## Elemental

    git clone https://gitlab.com/bootstrapcollaboration/elemental.git
    cd elemental
    mkdir build
    cd build
    cmake .. -DCMAKE_INSTALL_PREFIX=$HOME/install -DCMAKE_CXX_COMPILER=mpicxx -DCMAKE_C_COMPILER=mpicc
    make && make install
    cd ../..

## FLINT

    git clone https://github.com/flintlib/flint.git
    cd flint
    ./bootstrap.sh
    ./configure CC=mpicc CXX=mpicxx --prefix=$HOME/install --host=broadwell
    make
    make check 
    make install

`--host=broadwell` makes sure that the library isn't compiled with `avx512` instructions, which will crash the older nodes existing e.g. in `pi_poland` partition. Note also that `Autoconf`, `Automake`, and `libtool` modules should be loaded.

## MPSolve

    git clone https://github.com/robol/MPSolve.git
    cd MPSolve
    ./autogen.sh
    CC=mpicc CXX=mpicxx ./configure --prefix=$HOME/install --disable-dependency-tracking --disable-examples --disable-ui --disable-graphical-debugger --disable-documentation
    make && make install
    cd ..

## sdpb

    git clone https://github.com/davidsd/sdpb.git
    cd sdpb
    CXX=mpicxx ./waf configure --elemental-dir=$HOME/install --flint-dir=$HOME/install --mpsolve-dir=$HOME/install --prefix=$HOME/install/sdpb/master
    ./waf # -j 1
    ./test/run_all_tests.sh
    ./waf install
    cd ..

Note that computationally intensive processes (compilation, tests etc.) should be run on compute nodes and not on login
nodes. See [HPC documentation](https://docs.ycrc.yale.edu/clusters-at-yale/job-scheduling/) for more details.
