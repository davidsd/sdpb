# System info
    Red Hat Enterprise Linux release 8.8 (Ootpa)
    Linux login1.grace.ycrc.yale.edu 4.18.0-477.36.1.el8_8.x86_64 #1 SMP Thu Nov 9 08:12:18 EST 2023 x86_64 x86_64 x86_64 GNU/Linux

### Useful links:

- [Yale HPC documentation](https://docs.ycrc.yale.edu/clusters-at-yale/)
- [Run Jobs with Slurm](https://docs.ycrc.yale.edu/clusters-at-yale/job-scheduling/)
- [Applications & Software ](https://docs.ycrc.yale.edu/clusters-at-yale/applications/)

# Load modules

For compiling and/or running SDPB, you have to load modules first:

    module load Autoconf/2.71-GCCcore-12.2.0 Automake/1.16.5-GCCcore-12.2.0 Boost/1.83.0-GCC-12.2.0 CMake/3.24.3-GCCcore-12.2.0 GMP/6.2.1-GCCcore-12.2.0 libarchive/3.6.1-GCCcore-12.2.0 libtool/2.4.7-GCCcore-12.2.0 MPFR/4.2.0-GCCcore-12.2.0 OpenBLAS/0.3.21-GCC-12.2.0 OpenMPI/4.1.4-GCC-12.2.0 RapidJSON/1.1.0-GCCcore-12.2.0

You may run `module -t list` to view loaded modules,
and `module purge` to unload all modules.

You may also add this command to your `~/.bashrc` file, so that modules will load automatically.

# Use existing SDPB installation

## Choose SDPB version

SDPB is installed in `/home/vd289/install/sdpb-<VERSION>` folder,
where `<VERSION>` denotes specific version.

You may list all available versions via

    ls /home/vd289/install | grep sdpb

Fo example, `sdpb-master` is built from the latest [master](https://github.com/davidsd/sdpb/tree/master) branch (
run `sdpb --version` to see commit hash, e.g. `SDPB 2.5.1-130-g88b1c9ae`),
and `sdpb-3.0.0` is a stable [3.0.0](https://github.com/davidsd/sdpb/releases/tag/3.0.0) release.

Examples below are for `sdpb-master`.
You may replace it with another version, e.g. `sdpb-3.0.0`.

## Run SDPB

    /home/vd289/install/sdpb-master/bin/sdpb --help

### Batch script example

    sbatch /home/vd289/install/sdpb-master/share/sdpb_example.sh

This command submits `sdpb_example.sh` to
the [queueing system](https://docs.ycrc.yale.edu/clusters-at-yale/job-scheduling/).

`sdpb_example.sh` loads modules and runs `pmp2sdp`+`sdpb` for a simple problem.
See script code and comments for more details.

Script output is written to the log file in the current directory, e.g.:
`./sdpb_example.sh.26306151.out`.
SDPB output files are written to the `./out/` folder in the current directory.

# Build SDPB from sources

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

## sdpb

    git clone https://github.com/davidsd/sdpb.git
    cd sdpb
    CXX=mpicxx python3 ./waf configure --prefix=$HOME/install/sdpb-master --elemental-dir=$HOME/install --flint-dir=$HOME/install
    python3 ./waf # -j 1
    ./test/run_all_tests.sh
    python3 ./waf install
    cd ..