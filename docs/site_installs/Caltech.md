# System info

    CentOS Linux release 7.9.2009 (Core)
    Linux login2.cm.cluster 3.10.0-1160.53.1.el7.x86_64 #1 SMP Fri Jan 14 13:59:45 UTC 2022 x86_64 x86_64 x86_64 GNU/Linux

### Useful links:

- [Caltech HPC documentation](https://hpc.sites.caltech.edu/documentation)
- [How to run MPI jobs](https://hpc.sites.caltech.edu/documentation/slurm-commands)
- [Software and Modules](https://hpc.sites.caltech.edu/documentation/software-and-modules)

# Use existing SDPB installation

## Load modules

    module use /home/vdommes/install/modules
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

    sbatch $SDPB_HOME/share/sdpb_example.sh

submits the script to the [queueing system](https://hpc.sites.caltech.edu/documentation/slurm-commands).

`sdpb_example.sh` loads modules and runs `pmp2sdp`+`sdpb` for a simple problem.
See script code and comments for more details.

Script output is written to the log file in the current directory, e.g.:
`./sdpb_example.sh.26306151.out`.
SDPB output files are written to the `./out/` folder in the current directory.

# Build SDPB from sources

## Load modules

For compiling and/or running your custom SDPB installation, you have to load modules first:

    module load boost/1.84.0-gcc-11.3.1-zauawkb openmpi/5.0.1-gcc-13.2.0-xpoh5uw mpfr/4.2.0-gcc-13.2.0-yy2fkq5 libarchive/3.7.1-gcc-13.2.0-bzt555v openblas/0.3.23-gcc-11.3.1-cu3huj2 

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
    ./configure CC=mpicc CXX=mpicxx --disable-static --with-gmp=`pkg-config --variable=prefix gmp` --with-mpfr=`pkg-config --variable=prefix mpfr`--prefix=$HOME/install 
    make
    make check 
    make install
    cd ..

## MPSolve

    git clone https://github.com/robol/MPSolve.git
    cd MPSolve
    ./autogen.sh
    CC=mpicc CXX=mpicxx ./configure --prefix=$HOME/install --disable-dependency-tracking --disable-examples --disable-ui --disable-graphical-debugger --disable-documentation
    make && make install
    cd ..

## RapidJSON

    git clone https://github.com/Tencent/rapidjson.git
    cp -r rapidjson/include $HOME/install

## sdpb

    CXX=mpicxx ./waf configure --elemental-dir=$HOME/install --gmpxx-dir=/central/software9/spack/opt/spack/linux-rhel9-x86_64/gcc-13.2.0/gmp-6.2.1-lcnhysea5isaes4r547jt5nw2qru5ab7 --libarchive-dir=$(pkg-config --variable=prefix libarchive) --mpfr-dir=$(pkg-config --variable=prefix mpfr) --rapidjson-dir=$HOME/install --flint-dir=$HOME/install --mpsolve-dir=$HOME/install --boost-dir=$BOOST_ROOT --prefix=$HOME/install/sdpb/master
    ./waf # -j 1
    ./test/run_all_tests.sh
    ./waf install
    cd ..

Note that computationally intensive processes (compilation, tests etc.) should be run on compute nodes and not on login
nodes. See [HPC documentation](https://hpc.sites.caltech.edu/documentation/slurm-commands) for more details.
