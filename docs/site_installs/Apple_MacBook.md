# Apple MacBook

These instructions have been tested for MacBook Air M2 and for MacBook Air 2017 (Intel i5).

The only difference is that for ARM processors you have to set `-arch arm64` compiler flag for sdpb (see below).

# Build SDPB from sources

## Homebrew packages

If you don't have Homebrew package manager, install it following official instructions from https://brew.sh/. Add brew
to path following instructions shown after installing.

Then you can install packages required for SDPB:

    brew install gmp mpfr boost rapidjson libarchive openblas cmake open-mpi

You can see installation directory and another information for a package (e.g. `boost`)
by calling `brew info <package>` (e.g. `brew info boost`).

## Python 2

If Python 2 is missing, install it using macOS installer from https://www.python.org/downloads/release/python-2718/

You may try to skip this step, but note that waf fails with `bunzip2: Data integrity error when decompressing` when
using default Python 3, at least on some laptops.

## Elemental

    git clone https://gitlab.com/bootstrapcollaboration/elemental.git
    cd elemental
    mkdir build
    cd build
    cmake .. -DCMAKE_INSTALL_PREFIX=$HOME/install -DCMAKE_CXX_COMPILER=mpicxx -DCMAKE_C_COMPILER=mpicc
    make && make install
    cd ../..

## sdpb

### Checkout
    git clone https://github.com/davidsd/sdpb.git
    cd sdpb 

### Configure

    ./waf configure --elemental-dir=$HOME/install --boost-dir=/opt/homebrew/opt/boost/ --gmpxx-dir=/opt/homebrew/Cellar/gmp/6.3.0/ --mpfr-dir=/opt/homebrew/Cellar/mpfr/4.2.1/ --rapidjson-dir=/opt/homebrew/Cellar/rapidjson/1.1.0 --libarchive-dir=/opt/homebrew/Cellar/libarchive/3.7.2/ --prefix=$HOME/install/sdpb-master

If waf fails to find some package, e.g. `boost`, check the installation directory by calling, e.g. `brew info boost` and
update `--boost-dir` argument above.

The above `./waf configure` command works or x86 processors (e.g. Intel i5), but may fail for or ARM processors (e.g.
M2) with linker warnings `found architecture 'arm64', required architecture 'x86_64` in `build/config.log`.
In that case, you should set `-arch arm64` flag explicitly:

    CXXFLAGS="${CXXFLAGS} -arch arm64" LDFLAGS="${LDFLAGS} -arch arm64" ./waf configure --elemental-dir=$HOME/install --boost-dir=/opt/homebrew/opt/boost/ --gmpxx-dir=/opt/homebrew/Cellar/gmp/6.3.0/ --mpfr-dir=/opt/homebrew/Cellar/mpfr/4.2.1/ --rapidjson-dir=/opt/homebrew/Cellar/rapidjson/1.1.0 --libarchive-dir=/opt/homebrew/Cellar/libarchive/3.7.2/ --prefix=$HOME/install/sdpb-master

### Compile and install

    ./waf                   # build binaries
    ./build/sdpb --help     # run SDPB binary
    ./test/run_all_tests.sh # run tests to check correctness
    ./waf install           # install sdpb to --prefix=$HOME/install/sdpb-master
