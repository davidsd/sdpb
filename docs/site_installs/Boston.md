# System info

    AlmaLinux release 8.8 (Sapphire Caracal)
    Linux scc2 4.18.0-477.15.1.el8_8.x86_64 #1 SMP Mon Jul 24 04:58:00 EDT 2023 x86_64 x86_64 x86_64 GNU/Linux

### Useful links

- [BU Shared Computing Cluster Documentation](https://www.bu.edu/tech/support/research/system-usage/)
- [Running MPI jobs](https://www.bu.edu/tech/support/research/system-usage/running-jobs/)

# Load modules

For compiling and/or running SDPB, you have to load modules first:

    module load python3 gcc/12.2.0 openmpi/4.1.5 cmake/3.22.2 gmplib/6.2.1 mpfr/4.0.1 boost/1.79.0_gcc5+ openblas/0.3.23

You may run `module -t list` to view loaded modules,
and `module purge` to unload all modules.

You may also add this command to your `~/.bashrc` file, so that modules will load automatically.

# Use existing SDPB installation

## Choose SDPB version

SDPB is installed in `/usr2/collab/vdommes/install/sdpb-<VERSION>` folder,
where `<VERSION>` denotes specific version.

You may list all available versions via

    ls /usr2/collab/vdommes/install | grep sdpb

Fo example, `sdpb-master` is built from the latest [master](https://github.com/davidsd/sdpb/tree/master) branch (
run `sdpb --version` to see commit hash, e.g. `SDPB 2.5.1-130-g88b1c9ae`),
and `sdpb-2.6.1` is a stable [2.6.1](https://github.com/davidsd/sdpb/releases/tag/2.6.1) release.

Examples below are for `sdpb-master`.
You may replace it with another version, e.g. `sdpb-2.6.1`.
In that case, please refer
to [2.6.1 documentation](https://github.com/davidsd/sdpb/blob/2.6.1/docs/site_installs/Boston.md).

## Run SDPB

    /usr2/collab/vdommes/install/sdpb-master/bin/sdpb --help

### Batch script example

    qsub /usr2/collab/vdommes/install/sdpb-master/share/sdpb_example.sh

This command submits `sdpb_example.sh` to
the [queueing system](https://www.bu.edu/tech/support/research/system-usage/running-jobs/).

`sdpb_example.sh` loads modules and runs `sdp2input`+`sdpb` for a simple problem.
See script code and comments for more details.

Script output is written to the log file in the current directory, e.g.:
`./run_sdpb_example.sh.o8388881`.
SDPB output files are written to the `./out/` folder in the current directory.

# Build SDPB from sources

Use `RPATH` instead of `RUNPATH` in `mpicxx` linker, to fix shared library loading in SDPB:

    export OMPI_LDFLAGS="$(mpicxx --showme:link) -Wl,--disable-new-dtags"

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
    wget http://www.libarchive.org/downloads/libarchive-3.7.1.tar.xz
    tar -xf libarchive-3.7.1.tar.xz
    cd libarchive-3.7.1
    ./configure --prefix=$HOME/install
    make && make install
    cd ..

## sdpb
    git clone https://github.com/davidsd/sdpb.git
    cd sdpb 
    ./waf configure --elemental-dir=$HOME/install --rapidjson-dir=$HOME/install --libarchive-dir=$HOME/install --prefix=$HOME/install/sdpb-master
    ./waf # -j 1
    ./test/run_all_tests.sh
    ./waf install
    cd ..
