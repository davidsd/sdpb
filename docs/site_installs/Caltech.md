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

# Install SDPB

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

Install scalar_blocks
=============

Trilinos
--------
    git clone --branch trilinos-release-12-12-branch https://github.com/trilinos/Trilinos.git
    mkdir build
    cd build
    cmake -DCMAKE_CXX_COMPILER=g++ -DCMAKE_C_COMPILER=gcc -DTrilinos_ENABLE_Sacado=ON -DTrilinos_ENABLE_Kokkos=OFF -DTrilinos_ENABLE_Teuchos=OFF -DCMAKE_INSTALL_PREFIX=$HOME/install ..
    make && make install
    cd ../..

scalar_blocks
-------------
    git clone https://gitlab.com/bootstrapcollaboration/scalar_blocks.git
    cd scalar_blocks
    ./waf configure --prefix=$HOME/install --trilinos-dir=$HOME/install --eigen-incdir=/software/eigen-b3f3d4950030/
    ./waf # maybe -j 1
    ./waf install
    cd ..

Install blocks_3d
=========

fmt
---
    wget https://github.com/fmtlib/fmt/releases/download/6.2.1/fmt-6.2.1.zip
    mkdir build
    cd build
    cmake -DCMAKE_CXX_COMPILER=g++ -DCMAKE_INSTALL_PREFIX=$HOME/install ..
    make && make install
    cd ../..
  
blocks_3d
---------
    git clone git@gitlab.com:bootstrapcollaboration/blocks_3d.git
    cd blocks_3d
    ./waf configure --prefix=$HOME/install --eigen-incdir=/software/eigen-b3f3d4950030/ --fmt-dir=$HOME/install --fmt-libdir=$HOME/install/lib64
    ./waf # maybe -j
    ./waf install
    cd ..

Note: If you are trying to build the "profiling" branch of `blocks_3d`, you will need to load the following modules before attempting the above:

    module purge
    module load gcc/9.2.0
    module load boost/1_76_0_gcc-9.2.0
    module load python3/3.7.0
    module load eigen/eigen

Batch scripts
-------------
    /home/wlandry/sdpb/runs/TTTT_small.sh
    /home/wlandry/scalar_blocks/runs/scalar_blocks.sh
