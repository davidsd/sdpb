# System info

    AlmaLinux release 8.8 (Sapphire Caracal)
    Linux scc2 4.18.0-477.15.1.el8_8.x86_64 #1 SMP Mon Jul 24 04:58:00 EDT 2023 x86_64 x86_64 x86_64 GNU/Linux

### Useful links

- [BU Shared Computing Cluster Documentation](https://www.bu.edu/tech/support/research/system-usage/)
- [Running MPI jobs](https://www.bu.edu/tech/support/research/system-usage/running-jobs/)

# Use existing SDPB installation

## Load modules

    module use /home/vdommes/install/modules
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
    module load sdpb/3.1.0

## Run SDPB

An example script can be found at `$SDPB_HOME/share/sdpb_example.sh`.
The command

    qsub $SDPB_HOME/share/sdpb_example.sh

submits the script to the [queueing system](https://www.bu.edu/tech/support/research/system-usage/running-jobs/).

`sdpb_example.sh` loads modules and runs `pmp2sdp`+`sdpb` for a simple problem.
See script code and comments for more details.

Script output is written to the log file in the current directory, e.g.:
`./sdpb_example.o26306151`.
SDPB output files are written to the `./out/` folder in the current directory.

# Build SDPB from sources

## Load modules

For compiling and/or running your custom SDPB installation, you have to load modules first:

    module load python3 gcc/12.2.0 openmpi/4.1.5 cmake/3.22.2 gmplib/6.2.1 mpfr/4.2.1 boost/1.83.0 openblas/0.3.23

## Set environment variables

Use `RPATH` instead of `RUNPATH` in `mpicxx` linker, to fix shared library loading in SDPB:

    export OMPI_LDFLAGS="$(mpicxx --showme:link) -Wl,--disable-new-dtags"

## Elemental

    git clone https://gitlab.com/bootstrapcollaboration/elemental.git
    cd elemental
    mkdir build
    cd build
    GMPDIR=$SCC_GMPLIB_DIR MPFRDIR=$SCC_MPFR_DIR cmake .. -DCMAKE_INSTALL_PREFIX=$HOME/install -DCMAKE_CXX_COMPILER=mpicxx -DCMAKE_C_COMPILER=mpicc
    make && make install
    cd ../..

## FLINT

NB: when compiling FLINT on a login node, SDPB fails at runtime with
`SIGILL #0 __funlockfile #1 n_is_probabprime #2 comb` error.

You should compile it on a compute node,
e.g. [loading interactive session](https://www.bu.edu/tech/support/research/system-usage/running-jobs/interactive-jobs/)
via `qrsh`.

    git clone https://github.com/flintlib/flint.git
    cd flint
    ./bootstrap.sh
    ./configure --prefix=$HOME/install --disable-static --disable-pthread --with-gmp=$SCC_GMPLIB_DIR --with-mpfr=$SCC_MPFR_DIR CC=mpicc CXX=mpicxx
    make #-j1
    make check
    make install
    cd ..

## RapidJSON

    git clone https://github.com/Tencent/rapidjson.git
    cp -r rapidjson/include $HOME/install

## libarchive

    wget http://www.libarchive.org/downloads/libarchive-3.7.1.tar.xz
    tar -xf libarchive-3.7.1.tar.xz
    cd libarchive-3.7.1
    ./configure --prefix=$HOME/install
    make && make install
    cd ..

## MPSolve

    git clone https://github.com/robol/MPSolve.git
    cd MPSolve
    ./autogen.sh
    ./configure --prefix=$HOME/install --disable-dependency-tracking --disable-examples --disable-ui --disable-graphical-debugger --disable-documentation
    make && make install
    cd ..

## sdpb

    git clone https://github.com/davidsd/sdpb.git
    cd sdpb
    ./waf configure --elemental-dir=$HOME/install --rapidjson-dir=$HOME/install --libarchive-dir=$HOME/install --gmpxx-dir=$SCC_GMPLIB_DIR --mpfr-dir=$SCC_MPFR_DIR --flint-dir=$HOME/install --cblas-dir=$SCC_OPENBLAS_DIR --mpsolve-dir=$HOME/install --prefix=$HOME/install/sdpb/master
    ./waf # -j 1
    ./test/run_all_tests.sh
    ./waf install
    cd ..
