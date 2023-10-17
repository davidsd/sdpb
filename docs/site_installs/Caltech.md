# System info

    CentOS Linux release 7.9.2009 (Core)
    Linux login2.cm.cluster 3.10.0-1160.53.1.el7.x86_64 #1 SMP Fri Jan 14 13:59:45 UTC 2022 x86_64 x86_64 x86_64 GNU/Linux

### Useful links:

- [Caltech HPC documentation](https://hpc.sites.caltech.edu/documentation)
- [How to run MPI jobs](https://hpc.sites.caltech.edu/documentation/slurm-commands)
- [Software and Modules](https://hpc.sites.caltech.edu/documentation/software-and-modules)

# Load modules

    module load cmake/3.25.1 gcc/9.2.0 openmpi/4.1.5 boost/1_81_0_openmpi-4.1.1_gcc-9.2.0 eigen/eigen mpfr/4.0.2

You may run `module -t list` to view loaded modules,
and `module purge` to unload all modules.

# Use existing SDPB installation

SDPB executables built from the latest `master` branch can be found in `/home/vdommes/install/sdpb-master/bin/` folder,
e.g.

    /home/vdommes/install/sdpb-master/bin/sdpb --help

Stable release versions are also available, e.g.:

    /home/vdommes/install/sdpb-2.6.0/bin/sdpb --help

NB: remember to load modules before using SDPB.

# Build SDPB from sources

TODO: Note that `/usr/include/gmp.h` is present on login node, but may be absent on compute nodes. In that case one
should compile on a login node.

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
    cd ../..

## sdpb

    ./waf configure --elemental-dir=$HOME/install --rapidjson-dir=$HOME/install --libarchive-dir=$HOME/install --mpfr-dir=/central/software/mpfr/4.0.2 --prefix=$HOME/install/sdpb-master
    ./waf # I needed to do './waf -j 1' (single threaded) to get it to compile without crashing
    ./test/run_all_tests.sh
    ./waf install
    cd ..
