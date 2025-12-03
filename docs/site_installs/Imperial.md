# System info

    Red Hat Enterprise Linux release 8.5 (Ootpa)
    Linux login-c 4.18.0-348.20.1.el8_5.x86_64 #1 SMP Tue Mar 8 12:56:54 EST 2022 x86_64 x86_64 x86_64 GNU/Linux

### Useful links for Imperial College London HPC:

- [HPC documentation](https://icl-rcs-user-guide.readthedocs.io/en/latest/hpc/)
- [MPI jobs](https://icl-rcs-user-guide.readthedocs.io/en/latest/hpc/queues/mpi-jobs/)
- [Applications](https://icl-rcs-user-guide.readthedocs.io/en/latest/hpc/applications/)

# Use existing SDPB installation

## Load modules

    module use /rds/general/user/vdommes/install/modules
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

submits the script to the [queueing system](https://icl-rcs-user-guide.readthedocs.io/en/latest/hpc/queues/mpi-jobs/).

`sdpb_example.sh` loads modules and runs `pmp2sdp`+`sdpb` for a simple problem.
See script code and comments for more details.

Script output is written to the log file in the current directory, e.g.:
`./sdpb_example.sh.o26306151`.
SDPB output files are written to the `./out/` folder in the current directory.

# Build SDPB from sources

## Load modules

For compiling and/or running your custom SDPB installation, you have to load modules first:

    module load tools/prod Autotools/20220317-GCCcore-12.3.0 Boost/1.82.0-GCC-12.3.0 CMake/3.26.3-GCCcore-12.3.0 GCC/12.3.0 GMP/6.2.1-GCCcore-12.3.0 libarchive/3.6.2-GCCcore-12.3.0 METIS/5.1.0-GCCcore-12.3.0 MPFR/4.2.0-GCCcore-12.3.0 OpenBLAS/0.3.23-GCC-12.3.0 OpenMPI/4.1.5-GCC-12.3.0 pkg-config/0.29.2-GCCcore-12.3.0 Python/3.11.3-GCCcore-12.3.0 RapidJSON/1.1.0-20230928-GCCcore-12.3.0

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
    ./configure CC=mpicc CXX=mpicxx --prefix=$HOME/install
    make && make check && make install
    cd ..

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
    ./waf
    ./test/run_all_tests.sh
    ./waf install
    cd ..
