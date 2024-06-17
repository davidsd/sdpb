# System info

    Red Hat Enterprise Linux release 8.5 (Ootpa)
    Linux login-c 4.18.0-348.20.1.el8_5.x86_64 #1 SMP Tue Mar 8 12:56:54 EST 2022 x86_64 x86_64 x86_64 GNU/Linux

### Useful links for Imperial College London HPC:

- [HPC documentation](https://icl-rcs-user-guide.readthedocs.io/en/latest/hpc/)
- [MPI jobs](https://icl-rcs-user-guide.readthedocs.io/en/latest/hpc/queues/mpi-jobs/)
- [Applications](https://icl-rcs-user-guide.readthedocs.io/en/latest/hpc/applications/)

# Load modules

For compiling and/or running SDPB, you have to load modules first:

    module load cmake/3.18.2 oneapi/mpi/2021.4.0 gcc/11.2.0 boost/1.78.0 metis/5.1.0

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
and `sdpb-3.0.0` is a stable [3.0.0](https://github.com/davidsd/sdpb/releases/tag/3.0.0) release.

Examples below are for `sdpb-master`.
You may replace it with another version, e.g. `sdpb-3.0.0`.
In that case, please refer
to [3.0.0 documentation](https://github.com/davidsd/sdpb/blob/3.0.0/docs/site_installs/Imperial.md).

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

## GMP
We have to build GMP manually: FLINT requires at least GMP 6.2.1, but only GMP 5.0.1 is available in `module av gmp`.

    curl -o gmp-6.3.0.tar.xz https://gmplib.org/download/gmp/gmp-6.3.0.tar.xz
    tar -xf gmp-6.3.0.tar.xz
    cd gmp-6.3.0
    ./configure --enable-cxx --prefix=$HOME/install
    make #-j1
    make check
    make install
    cd ..

You may additionally tune up performance following instructions at https://gmplib.org/manual/Performance-optimization

## MPFR
We have to build GMP manually: FLINT requires at least MPFR 4.1.0, but only MPFR 3.1.1 is available in `module av mpfr`.

    curl -o mpfr-4.2.1.tar.xz https://www.mpfr.org/mpfr-current/mpfr-4.2.1.tar.xz
    tar -xf mpfr-4.2.1.tar.xz
    cd mpfr-4.2.1
    ./configure --with-gmp=$HOME/install --prefix=$HOME/install
    make #-j1
    make check
    make install
    cd ..

## Elemental

    git clone https://gitlab.com/bootstrapcollaboration/elemental.git
    cd elemental
    mkdir build
    cd build    
    cmake .. -DCMAKE_INSTALL_PREFIX=$HOME/install -DCMAKE_CXX_COMPILER=mpicxx -DCMAKE_C_COMPILER=mpicc -DGMP_INCLUDES=$HOME/install/include -DGMP_LIBRARIES=$HOME/install/lib/libgmp.so -DGMPXX_INCLUDES=$HOME/install/include -DGMPXX_LIBRARIES=$HOME/install/lib/libgmpxx.so -DMPFR_INCLUDES=$HOME/install/include -DMPFR_LIBRARIES=$HOME/install/lib/libmpfr.so
    make && make install
    cd ../..

## FLINT

    git clone https://github.com/flintlib/flint.git
    cd flint
    ./bootstrap.sh
    ./configure --disable-static --prefix=$HOME/install --with-gmp=$HOME/install --with-mpfr=$HOME/install CC=mpicc CXX=mpicxx
    make #-j1
    make check
    make install 
    cd ..

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
    python3 ./waf configure --gmpxx-dir=$HOME/install --mpfr-dir=$HOME/install --elemental-dir=$HOME/install --rapidjson-dir=$HOME/install --libarchive-dir=$HOME/install --flint-dir=$HOME/install --prefix=$HOME/install/sdpb-master
    python3 ./waf # -j 1
    ./test/run_all_tests.sh mpirun
    python3 ./waf install
    cd ..
