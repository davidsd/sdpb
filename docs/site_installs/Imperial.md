# System info

    Red Hat Enterprise Linux release 8.5 (Ootpa)
    Linux login-c 4.18.0-348.20.1.el8_5.x86_64 #1 SMP Tue Mar 8 12:56:54 EST 2022 x86_64 x86_64 x86_64 GNU/Linux

### Useful links for Imperial College London HPC:

- [HPC documentation](https://wiki.imperial.ac.uk/display/HPC)
- [MPI jobs](https://wiki.imperial.ac.uk/display/HPC/MPI+Jobs)
- [Applications](https://wiki.imperial.ac.uk/display/HPC/Applications)

# Load modules

For compiling and/or running SDPB, you have to load modules first:

    module load cmake/3.18.2 oneapi/mpi/2021.4.0 gcc/11.2.0 boost/1.78.0 gmp/5.0.1 mpfr/3.1.1 metis/5.1.0

You may run `module -t list` to view loaded modules,
and `module purge` to unload all modules.

You may also add this command to your `~/.bashrc` file, so that modules will load automatically.

# Use existing SDPB installation

## Choose SDPB version

SDPB is installed in `/rds/general/user/vdommes/home/install/sdpb-<VERSION>` folder,
where `<VERSION>` denotes specific version.

You may list all available versions via

    ls /rds/general/user/vdommes/home/install | grep sdpb

Fo example, `sdpb-master` is built from the latest [master](https://github.com/davidsd/sdpb/tree/master) branch (
run `sdpb --version` to see commit hash, e.g. `SDPB 2.5.1-130-g88b1c9ae`),
and `sdpb-2.6.1` is a stable [2.6.1](https://github.com/davidsd/sdpb/releases/tag/2.6.1) release.

Examples below are for `sdpb-master`.
You may replace it with another version, e.g. `sdpb-2.6.1`.
In that case, please refer
to [2.6.1 documentation](https://github.com/davidsd/sdpb/blob/2.6.1/docs/site_installs/Imperial.md).

## Run SDPB

    /rds/general/user/vdommes/home/install/sdpb-master/bin/sdpb --help

### Batch script example

    qsub /rds/general/user/vdommes/home/install/sdpb-master/share/sdpb_example.sh

This command submits `sdpb_example.sh` to
the [queueing system](https://wiki.imperial.ac.uk/display/HPC/Queueing+System).

`sdpb_example.sh` loads modules and runs `pmp2sdp`+`sdpb` for a simple problem.
See script code and comments for more details.

Script output is written to the log file in the current directory, e.g.:
`./run_sdpb_example.sh.o8388881`.
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

## RapidJSON

    git clone https://github.com/Tencent/rapidjson.git
    cp -r rapidjson/include $HOME/install

## libarchive

    git clone -b v3.6.2 https://github.com/libarchive/libarchive.git
    cd libarchive
    ./build/autogen.sh
    ./configure --prefix=$HOME/install
    make && make install
    cd ../..

## sdpb

    git clone https://github.com/davidsd/sdpb.git
    cd sdpb 
    python3 ./waf configure --mpfr-dir=$MPFR_HOME --elemental-dir=$HOME/install --rapidjson-dir=$HOME/install --libarchive-dir=$HOME/install --prefix=$HOME/install/sdpb-master
    python3 ./waf # -j 1
    ./test/run_all_tests.sh mpirun
    python3 ./waf install
    cd ..
