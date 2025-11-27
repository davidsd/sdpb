# System info

    Rocky Linux release 8.8 (Green Obsidian)
    Linux exp-1-33 4.18.0-477.15.1.el8_8.x86_64 #1 SMP Wed Jun 28 15:04:18 UTC 2023 x86_64 x86_64 x86_64 GNU/Linux

### Useful links

- [Expanse HPC User Guide](https://www.sdsc.edu/support/user_guides/expanse.html)
- [Running Jobs on Expanse](https://www.sdsc.edu/systems/expanse/user_guide.html#narrow-wysiwyg-7)

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
    module load sdpb/3.0.0

## Run SDPB

An example script can be found at `$SDPB_HOME/share/sdpb_example.sh`.
The command

    sbatch $SDPB_HOME/share/sdpb_example.sh

submits the script to the [queueing system](https://www.sdsc.edu/systems/expanse/user_guide.html#narrow-wysiwyg-7).

`sdpb_example.sh` loads modules and runs `pmp2sdp`+`sdpb` for a simple problem.
See script code and comments for more details.

Script output is written to the log file in the current directory, e.g.:
`./sdpb_example.sh.26306151.out`.
SDPB output files are written to the `./out/` folder in the current directory.

# Build SDPB from sources

## Load modules

For compiling and/or running your custom SDPB installation, you have to load modules first:

    module load cpu/0.15.4 gcc/10.2.0 openmpi/4.0.4 gmp/6.1.2 mpfr/4.0.2 cmake/3.18.2 openblas/dynamic/0.3.7 boost/1.74.0 slurm

## Elemental

    git clone https://gitlab.com/bootstrapcollaboration/elemental.git
    cd elemental
    mkdir build
    cd build
    cmake .. -DCMAKE_CXX_COMPILER=mpicxx -DCMAKE_C_COMPILER=mpicc -DCMAKE_Fortran_COMPILER=mpif90 -DCMAKE_INSTALL_PREFIX=$HOME/install -DMATH_LIBS="-L$BLASDIR -lopenblas"
    make -j 64 && make install
    cd ../..

## FLINT

    git clone https://github.com/flintlib/flint.git -b v3.0.1
    cd flint
    ./bootstrap.sh
    ./configure CC=mpicc CXX=mpicxx --disable-static --disable-gmp-internals --prefix=$HOME/install 
    make
    make check 
    make install

Note that you have to checkout `v3.0.1` because starting from v3.1.0 FLINT requires GMP 6.2.1 or later.
The flag `--disable-gmp-internals` is required to prevent `mpn_gcd_11 not found` error.

## libarchive

    wget http://www.libarchive.org/downloads/libarchive-3.7.1.tar.xz
    tar -xf libarchive-3.7.1.tar.xz
    cd libarchive-3.7.1
    ./configure --prefix=$HOME/install
    make -j 64 && make install
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

    git clone https://github.com/davidsd/sdpb.git
    cd sdpb
    CXX=mpicxx ./waf configure --elemental-dir=$HOME/install --flint-dir=$HOME/install --rapidjson-dir=$HOME/install --boost-dir=$HOME/install  --libarchive-dir=$HOME/install --mpsolve-dir=$HOME/install --prefix=$HOME/install/sdpb/master
    ./waf # -j 1
    ./test/run_all_tests.sh
    ./waf install
    cd ..

Note that computationally intensive processes (compilation, tests etc.) should be run on compute nodes and not on login
nodes. See [Expanse User Guide](https://www.sdsc.edu/systems/expanse/user_guide.html#narrow-wysiwyg-7) for more details.
