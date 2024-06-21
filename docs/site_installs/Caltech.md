# System info

    CentOS Linux release 7.9.2009 (Core)
    Linux login2.cm.cluster 3.10.0-1160.53.1.el7.x86_64 #1 SMP Fri Jan 14 13:59:45 UTC 2022 x86_64 x86_64 x86_64 GNU/Linux

### Useful links:

- [Caltech HPC documentation](https://hpc.sites.caltech.edu/documentation)
- [How to run MPI jobs](https://hpc.sites.caltech.edu/documentation/slurm-commands)
- [Software and Modules](https://hpc.sites.caltech.edu/documentation/software-and-modules)

# Load modules

For compiling and/or running SDPB, you have to load modules first:

    module load boost/1.84.0-gcc-11.3.1-zauawkb gcc/13.2.0-gcc-13.2.0-w55nxkl openmpi/5.0.1-gcc-13.2.0-xpoh5uw libarchive/3.7.1-gcc-13.2.0-bzt555v openblas/0.3.23-gcc-13.2.0-3csynno

You may run `module -t list` to view loaded modules,
and `module purge` to unload all modules.

You may also add this command to your `~/.bashrc` file, so that modules will load automatically.

# Use existing SDPB installation

## Choose SDPB version

SDPB is installed in `/central/home/vdommes/install/sdpb-<VERSION>` folder,
where `<VERSION>` denotes specific version.

You may list all available versions via

    ls /central/home/vdommes/install | grep sdpb

Fo example, `sdpb-master` is built from the latest [master](https://github.com/davidsd/sdpb/tree/master) branch (
run `sdpb --version` to see commit hash, e.g. `SDPB 2.5.1-130-g88b1c9ae`),
and `sdpb-3.0.0` is a stable [3.0.0](https://github.com/davidsd/sdpb/releases/tag/3.0.0) release.

Examples below are for `sdpb-master`.
You may replace it with another version, e.g. `sdpb-3.0.0`.
In that case, please refer
to [3.0.0 documentation](https://github.com/davidsd/sdpb/blob/3.0.0/docs/site_installs/Caltech.md).

## Run SDPB

    /central/home/vdommes/install/sdpb-master/bin/sdpb --help

### Batch script example

    sbatch /central/home/vdommes/install/sdpb-master/share/sdpb_example.sh

This command submits `sdpb_example.sh` to
the [queueing system](https://hpc.sites.caltech.edu/documentation/slurm-commands).

`sdpb_example.sh` loads modules and runs `pmp2sdp`+`sdpb` for a simple problem.
See script code and comments for more details.

Script output is written to the log file in the current directory, e.g.:
`./sdpb_example.sh.26306151.out`.
SDPB output files are written to the `./out/` folder in the current directory.

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

## FLINT

    git clone https://github.com/flintlib/flint.git
    cd flint
    ./bootstrap.sh
    ./configure CC=mpicc CXX=mpicxx --disable-static --with-gmp=`pkg-config --variable=prefix gmp` --with-mpfr=`pkg-config --variable=prefix mpfr`--prefix=$HOME/install 
    make
    make check 
    make install

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
    CXX=mpicxx ./waf configure --libarchive-dir=`pkg-config --variable=prefix libarchive` --elemental-dir=$HOME/install --rapidjson-dir=$HOME/install --flint-dir=$HOME/install 
    --prefix=$HOME/install/sdpb-master
    ./waf # -j 1
    ./test/run_all_tests.sh
    ./waf install
    cd ..
