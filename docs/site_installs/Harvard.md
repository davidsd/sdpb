# System info

    Rocky Linux release 8.7 (Green Obsidian)
    Linux holylogin02.rc.fas.harvard.edu 4.18.0-425.10.1.el8_7.x86_64 #1 SMP Thu Jan 12 16:32:13 UTC 2023 x86_64 x86_64 x86_64 GNU/Linux

### Useful links

- [Harvard FAS RC User Guide](https://docs.rc.fas.harvard.edu/)
- [Running MPI jobs](https://docs.rc.fas.harvard.edu/kb/running-jobs/)

# Load modules

For compiling and/or running SDPB, you have to load modules first:

    module load python gcc/12.2.0-fasrc01 cmake/3.25.2-fasrc01 openmpi/4.1.5-fasrc01 gmp/6.2.1-fasrc01 mpfr/4.2.0-fasrc01

You may run `module -t list` to view loaded modules,
and `module purge` to unload all modules.

You may also add this command to your `~/.bashrc` file, so that modules will load automatically.

# Use existing SDPB installation

## Choose SDPB version

SDPB is installed in `/n/home02/vdommes/install/sdpb-<VERSION>` folder,
where `<VERSION>` denotes specific version.

You may list all available versions via

    ls /n/home02/vdommes/install | grep sdpb

Fo example, `sdpb-master` is built from the latest [master](https://github.com/davidsd/sdpb/tree/master) branch (
run `sdpb --version` to see commit hash, e.g. `SDPB 2.5.1-130-g88b1c9ae`),
and `sdpb-2.7.0` is a stable [2.7.0](https://github.com/davidsd/sdpb/releases/tag/2.7.0) release.

Examples below are for `sdpb-master`.
You may replace it with another version, e.g. `sdpb-2.7.0`.
In that case, please refer
to [2.7.0 documentation](https://github.com/davidsd/sdpb/blob/2.7.0/docs/site_installs/Harvard.md).

## Run SDPB

    /n/home02/vdommes/install/sdpb-master/bin/sdpb --help

### Batch script example

    sbatch /n/home02/vdommes/install/sdpb-master/share/sdpb_example.sh

This command submits `sdpb_example.sh` to
the [queueing system](https://docs.rc.fas.harvard.edu/kb/running-jobs/).

`sdpb_example.sh` loads modules and runs `pmp2sdp`+`sdpb` for a simple problem.
See script code and comments for more details.

Script output is written to the log file in the current directory, e.g.:
`./sdpb_example.sh.26306151.out`.
SDPB output files are written to the `./out/` folder in the current directory.

# Build SDPB from sources

Use `RPATH` instead of `RUNPATH` in `mpicxx` linker, to fix shared library loading in SDPB:

    export OMPI_LDFLAGS="$(mpicxx --showme:link) -Wl,--disable-new-dtags"

## Boost

    wget https://boostorg.jfrog.io/artifactory/main/release/1.83.0/source/boost_1_83_0.tar.bz2
    tar -jxf boost_1_83_0.tar.bz2
    cd boost_1_83_0
    ./bootstrap.sh --prefix=$HOME/install --without-libraries=python
    ./b2 -j 16 --prefix=$HOME/install
    ./b2 --prefix=$HOME/install install
    cd ..

## Elemental

    git clone https://gitlab.com/bootstrapcollaboration/elemental.git
    mkdir build
    cd build
    cmake .. -DCMAKE_INSTALL_PREFIX=$HOME/install -DCMAKE_CXX_COMPILER=mpicxx -DCMAKE_C_COMPILER=mpicc
    make && make install

## RapidJSON

    git clone https://github.com/Tencent/rapidjson.git
    cp -r rapidjson/include $HOME/install

## libarchive

    wget http://www.libarchive.org/downloads/libarchive-3.7.1.tar.xz
    tar -xf libarchive-3.7.1.tar.xz
    cd libarchive-3.7.1
    ./configure --prefix=$HOME/install
    make -j 16 && make install
    cd ..

## sdpb

    git clone https://github.com/davidsd/sdpb.git
    cd sdpb
    ./waf configure --boost-dir=$HOME/install --elemental-dir=$HOME/install --rapidjson-dir=$HOME/install --libarchive-dir=$HOME/install --prefix=$HOME/install/sdpb-master
    ./waf # -j 1
    ./test/run_all_tests.sh
    ./waf install
