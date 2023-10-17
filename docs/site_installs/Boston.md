# System info

    AlmaLinux release 8.8 (Sapphire Caracal)
    Linux scc2 4.18.0-477.15.1.el8_8.x86_64 #1 SMP Mon Jul 24 04:58:00 EDT 2023 x86_64 x86_64 x86_64 GNU/Linux

### Useful links

- [BU Shared Computing Cluster Documentation](https://www.bu.edu/tech/support/research/system-usage/)
- [Running MPI jobs](https://www.bu.edu/tech/support/research/system-usage/running-jobs/)

# Load modules

    module load python3 gcc/12.2.0 openmpi/4.1.5 cmake/3.22.2 gmplib/6.2.1 mpfr/4.0.1 boost/1.79.0_gcc5+ openblas/0.3.23

You may run `module -t list` to view loaded modules,
and `module purge` to unload all modules.

# Use existing SDPB installation

SDPB executables built from the latest `master` branch can be found in `/usr2/collab/vdommes/install/sdpb-master/bin/`
folder, e.g.

    /usr2/collab/vdommes/install/sdpb-master/bin/sdpb --help

Stable release versions are also available, e.g.:

    /usr2/collab/vdommes/install/sdpb-2.6.0/bin/sdpb --help

NB: remember to load modules before using SDPB.

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
    ./waf # I needed to do './waf -j 1' (single threaded) to get it to compile without crashing
    ./test/run_all_tests.sh
    ./waf install
    cd ..
