This guide is for building SDPB.  To just run SDPB, it may be easier
to use [Docker](docs/Docker.md) or [Singularity](docs/Singularity.md).

* [Requirements](#requirements)
* [Installation](#installation)

# Requirements

SDPB requires

- A modern C++ compiler with C++ 14 support.  SDPB has been tested with
  GCC 5.2.0, 7.1.0, and 8.2.0.

- An MPI implementation such as [OpenMPI](https://www.open-mpi.org/) or [MPICH](https://www.mpich.org/)

- [Boost C++ Libraries](http://www.boost.org/) Please be sure that the
  boost library is compiled with the same C++ compiler you are using.

- [The GNU Multiprecision Library](https://gmplib.org/).  Be sure to
  enable C++ support by configuring with the option `--enable-cxx`.

- [The GNU MPFR Library](https://www.mpfr.org/)

- [libxml2](http://www.xmlsoft.org/).  Only the C library is required.

- [RapidJSON](http://rapidjson.org/)

- A BLAS library such as [OpenBLAS](http://www.openblas.net/)

- [Python](https://python.org) (2.7 or later).

- [CMake](https://cmake.org/)

SDPB has only been tested on Linux (Debian buster and Centos 7).  On
Centos 7, the system compiler (gcc 4.8.5) is too old to support
C++ 14.  So you will have to install a newer compiler and Boost.  The
system versions of GMP, MPFR, and libxml2 have been tested to work.

In principle, SDPB should be installable on Mac OS X using a package
manager such as [Homebrew](https://brew.sh).

# Installation

1. Download the the fork of [Elemental](https://gitlab.com/bootstrapcollaboration/elemental)

        git clone https://gitlab.com/bootstrapcollaboration/elemental.git

2. Make the build directory and cd into it.

        mkdir -p elemental/build
        cd elemental/build

3. Configure Elemental.  This can be rather tricky.  You may have to specify where the Boost, GMP, and BLAS libraries are.  If you are using [modules](http://modules.sourceforge.net/), you may have to load modules.  If these define appropriate environment variables, then the configure should be relatively easy.  For example, on Yale's Grace cluster, commands that work are
   
        module load Langs/GCC/5.2.0 Libs/Boost Libs/OpenBLAS MPI/OpenMPI/1.8.1-gcc Tools/Cmake
        export CXX=mpicxx
        export CC=mpicc
        cmake .. -DCMAKE_INSTALL_PREFIX=$HOME/project/install

    On Harvard's Odyssey3 cluster, the default Boost module does not
    work with the most up-to-date compiler.  Using a different
    compiler with a different Boost module does work.
    
        module load cmake gcc/7.1.0-fasrc01 openmpi OpenBLAS boost/1.63.0-fasrc02 libxml2/2.7.8-fasrc03
        export CXX=mpicxx
        export CC=mpicc
        cmake .. -DCMAKE_INSTALL_PREFIX=$HOME/install

    On a Debian Buster laptop, use sudo or the root account to
    install these packages
    
        apt-get install openmpi-bin libopenmpi-dev libgmp-dev libmpfr-dev libmpfrc++-dev libboost-all-dev g++ cmake libopenblas-dev libxml2-dev git libmetis-dev pkg-config rapidjson-dev
    
    and then configure as a normal user
    
        cmake .. -DCMAKE_CXX_COMPILER=mpicxx -DCMAKE_C_COMPILER=mpicc -DCMAKE_INSTALL_PREFIX=/home/boo/qft/src/elemental/install

4. Build and install Elemental

        make && make install

5. Download the Elemental version of SDPB

        cd ../..
        git clone https://github.com/davidsd/sdpb
        cd sdpb

6. Configure the project using the included version of [waf](https://waf.io).  Waf can usually find libraries that are in system directories or have corresponding environment variables.  You may have to explicitly tell it where other libraries are.  If you are having problems, running `./waf --help` will give you a list of options.
   
   On Yale, a working command is

        ./waf configure --elemental-dir=$HOME/project/install

    On Harvard, a working command is

        ./waf configure --elemental-dir=$HOME/install

    and on Debian buster, it is

        ./waf configure --elemental-dir=$HOME/install
    
7. Type `./waf` to build the executable in `build/sdpb`.  This will create four executables in the `build/` directory: `pvm2sdp`, `sdp2blocks`, `block_grid_mapping`, and `sdpb`. Running
   
        ./build/sdpb --help

should give you a usage message.  HPC systems may require you to run
all jobs through a batch system.  You will need to consult the
documentation for your individual system.
