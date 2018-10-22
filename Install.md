
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

- [The GNU MPC Library](http://www.multiprecision.org/mpc)

- [libxml2](http://www.xmlsoft.org/).  Only the C library is required,
  so this should work with the version that comes with the operating
  system.

- A BLAS library such as [OpenBLAS](http://www.openblas.net/)

- [Python](https://python.org) (2.7 or later).

- [CMake](https://cmake.org/)

SDPB has only been tested on Linux (Debian buster and Centos 7).  In
principle, it should be installable on Mac OS X using a package
manager such as [Homebrew](https://brew.sh).


# Installation

1. Download the `no_warnings` branch of the fork of [Elemental](https://gitlab.com/bootstrapcollaboration/elemental)

    `git clone --branch no_warnings https://gitlab.com/bootstrapcollaboration/elemental`
    
2. Make the build directory and cd into it.

    `mkdir -p elemental/build`
    `cd elemental/build`
    
3. Configure Elemental.  This can be rather tricky.  This requires
   specifying where the Boost, GMP, MPFR, MPC, and the BLAS libraries
   are.  If you are using [modules](http://modules.sourceforge.net/),
   you may have to load modules.  For example, on Yale's Grace
   cluster, the commands that work are
   
   `module load Langs/GCC/5.2.0 Libs/Boost Libs/OpenBLAS MPI/OpenMPI/1.8.1-gcc Tools/Cmake Libs/GMP Libs/MPFR Libs/MPC`
   `export CXX=mpicxx`
   `export CC=mpicc`
   `cmake .. -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=$HOME/project/elemental/libelemental_bin -DGMP_INCLUDES=$GMP_INCLUDE -DGMP_LIBRARIES=$GMP_LIB/libgmp.so -DMPFR_INCLUDES=$MPFR_INCLUDE -DMPFR_LIBRARIES=$MPFR_LIB/libmpfr.so -DMPC_INCLUDES=$MPC_INCLUDE -DMPC_LIBRARIES=$MPC_LIB/libmpc.so -DEL_DISABLE_QD=ON -DEL_DISABLE_PARMETIS=true -DGMPXX_INCLUDES=$GMP_INCLUDE -DGMPXX_LIBRARIES=$GMP_LIB/libgmpxx.so -DOpenBLAS=$OPENBLAS_LIB/libopenblas.so -DGFORTRAN_LIB=/gpfs/apps/hpc/Langs/GCC/5.2.0/lib64/libgfortran.so`

    On Harvard's Odyssey3 cluster, the commands that work are
    
    `module load cmake gcc/7.1.0-fasrc01 openmpi OpenBLAS boost libxml2/2.7.8-fasrc03`
    `export CXX=mpicxx`
    `export CC=mpicc`
    `cmake .. -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=$HOME/elemental/libelemental_bin -DGMP_INCLUDES=$GMP_INCLUDE -DGMP_LIBRARIES=$GMP_LIB/libgmp.so -DMPFR_INCLUDES=$MPFR_INCLUDE -DMPFR_LIBRARIES=$MPFR_LIB/libmpfr.so -DMPC_INCLUDES=$MPC_INCLUDE -DMPC_LIBRARIES=$MPC_LIB/libmpc.so -DEL_DISABLE_QD=ON -DEL_DISABLE_PARMETIS=true -DGMPXX_INCLUDES=$GMP_INCLUDE -DGMPXX_LIBRARIES=$GMP_LIB/libgmpxx.so -DOpenBLAS=$OPENBLAS_LIB/libopenblas.so -DGFORTRAN_LIB=/n/helmod/apps/centos7/Core/gcc/7.1.0-fasrc01/lib64/libgfortran.so`

    On a Debian Buster laptop, use sudo or the root account to
    install these packages
    
    `apt-get install cmake, libopenblas-dev, openmpi-bin, libgmp-dev, libmpc-dev, libmpfr-dev, libboost-dev, libxml2-dev`
    
    and then configure as a normal user
    
    `export CXX=mpicxx`
    `export CC=mpicc`
    `cmake .. -DCMAKE_BUILD_TYPE=Release -DCMAKE_Fortran_COMPILER=mpif90 -DCMAKE_INSTALL_PREFIX=$HOME/elemental/libelemental_bin -DGFORTRAN_LIB=/usr/lib/gcc/x86_64-linux-gnu/8/libgfortran.so -DEL_DISABLE_QD=ON -DEL_DISABLE_PARMETIS=true`

4. Build and install Elemental

    `make && make install`

5. Download the Elemental version of SDPB

    `cd ../..`
    `git clone --branch elemental https://wlandry:chie6ieM@github.com/davidsd/sdpb`
    `cd sdpb`

6. Configure the project using the included version of
   [waf](https://waf.io).  Waf can usually find libraries that are in
   system directories, but it needs direction for everything else is.
   If you are having problems, running `python ./waf --help` will give
   you a list of options.
   
   On Yale, a working command is

    `python ./waf configure --elemental-incdir=$HOME/elemental/libelemental_bin/include --elemental-libdir="$HOME/project/elemental/libelemental_bin/lib64 $HOME/project/elemental/libelemental_bin/lib" --boost-includes=$BOOST_INCLUDE --boost-libs=$BOOST_LIB`

    On Harvard, a working command is

    `./waf configure --gmpxx-incdir=$GMP_INCLUDE --gmpxx-libdir=$GMP_LIB --elemental-libdir="$HOME/elemental/libelemental_bin/lib64 $HOME/elemental/libelemental_bin/lib" --elemental-incdir=$HOME/elemental/libelemental_bin/include --boost-includes=$BOOST_INCLUDE --boost-libs=$BOOST_LIB --libxml2-libdir=$LIBXML2_LIB --libxml2-incdir=$LIBXML2_INCLUDE`

    and on Debian buster, it is

    `python ./waf configure --elemental-dir=$HOME/elemental/libelemental_bin`
    

3. Type `python ./waf` to build the executable in `build/sdpb`.

