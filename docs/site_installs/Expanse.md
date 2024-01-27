# System info

    Rocky Linux release 8.8 (Green Obsidian)
    Linux exp-1-33 4.18.0-477.15.1.el8_8.x86_64 #1 SMP Wed Jun 28 15:04:18 UTC 2023 x86_64 x86_64 x86_64 GNU/Linux

### Useful links

- [Expanse HPC User Guide](https://www.sdsc.edu/support/user_guides/expanse.html)
- [Running MPI jobs](https://www.sdsc.edu/support/user_guides/expanse.html#running)

# Load modules

For compiling and/or running SDPB, you have to load modules first:

    module load cpu/0.15.4 gcc/10.2.0 openmpi/4.0.4 gmp/6.1.2 mpfr/4.0.2 cmake/3.18.2 openblas/dynamic/0.3.7 boost/1.74.0 slurm

You may run `module -t list` to view loaded modules,
and `module purge` to unload all modules.

You may also add this command to your `~/.bashrc` file, so that modules will load automatically.

# Use existing SDPB installation

## Choose SDPB version

SDPB is installed in `/home/vdommes/install/sdpb-<VERSION>` folder,
where `<VERSION>` denotes specific version.

You may list all available versions via

    ls /home/vdommes/install | grep sdpb

Fo example, `sdpb-master` is built from the latest [master](https://github.com/davidsd/sdpb/tree/master) branch (
run `sdpb --version` to see commit hash, e.g. `SDPB 2.5.1-130-g88b1c9ae`),
and `sdpb-2.6.1` is a stable [2.6.1](https://github.com/davidsd/sdpb/releases/tag/2.6.1) release.

Examples below are for `sdpb-master`.
You may replace it with another version, e.g. `sdpb-2.6.1`.
In that case, please refer
to [2.6.1 documentation](https://github.com/davidsd/sdpb/blob/2.6.1/docs/site_installs/Expanse.md).

## Run SDPB

    /home/vdommes/install/sdpb-master/bin/sdpb --help

### Batch script example

    sbatch /home/vdommes/install/sdpb-master/share/sdpb_example.sh

This command submits `sdpb_example.sh` to
the [queueing system](https://www.sdsc.edu/support/user_guides/expanse.html#running).

`sdpb_example.sh` loads modules and runs `pmp2sdp`+`sdpb` for a simple problem.
See script code and comments for more details.

Script output is written to the log file in the current directory, e.g.:
`./sdpb_example.sh.26306151.out`.
SDPB output files are written to the `./out/` folder in the current directory.

# Build SDPB from sources

TODO: on compute nodes compiler fails to find some libs

## Elemental

    git clone https://gitlab.com/bootstrapcollaboration/elemental.git
    cd elemental
    mkdir build
    cd build
    cmake .. -DCMAKE_CXX_COMPILER=mpicxx -DCMAKE_C_COMPILER=mpicc -DCMAKE_Fortran_COMPILER=mpif90 -DCMAKE_INSTALL_PREFIX=$HOME/install -DMATH_LIBS="-L$BLASDIR -lopenblas"
    make -j 64 && make install
    cd ../..

## libarchive

    wget http://www.libarchive.org/downloads/libarchive-3.7.1.tar.xz
    tar -xf libarchive-3.7.1.tar.xz
    cd libarchive-3.7.1
    ./configure --prefix=$HOME/install
    make -j 64 && make install
    cd ..

## RapidJSON

    git clone https://github.com/Tencent/rapidjson.git
    cp -r rapidjson/include $HOME/install

## sdpb

    git clone https://github.com/davidsd/sdpb.git
    cd sdpb
    CXX=mpicxx ./waf configure --prefix=$HOME/install/sdpb-master --elemental-dir=$HOME/install --rapidjson-dir=$HOME/install --boost-dir=$HOME/install  --libarchive-dir=$HOME/install
    ./waf # -j 1
    ./test/run_all_tests.sh
    ./waf install
    cd ..
