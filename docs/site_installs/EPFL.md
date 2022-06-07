# Helvetios site install

## Site packages

```bash
module load gcc/8.4.0 mvapich2/2.3.4 cmake/3.16.5
```

## Local directories

For this setup:

```bash
export SDPBDIR=/work/fsl/sdpb-local
mkdir -p $SDPBDIR
# chmod 775 $SDPBDIR
```

### `env`

```bash
export CC=mpicc
export CXX=mpic++
export CMAKE_C_COMPILER=mpicc
export CMAKE_CXX_COMPILER=mpic++
export CMAKE_INSTALL_PREFIX=$SDPBDIR
export CMAKE_PREFIX_PATH=$SDPBDIR:$CMAKE_PREFIX_PATH
```

## Upstream packages

### Boost

Version [1.79.0](https://boostorg.jfrog.io/artifactory/main/release/1.79.0/source/boost_1_79_0.tar.gz) used.

That one's bit hacky (PR fixing an issue on SDPB site will be done).

```bash
./bootstrap.sh --prefix=$SDPBDIR --without-libraries=python
echo using mpi \; >> project-config.jam
./b2 install --prefix=$SDPBDIR
```

And "a little patch":

```bash
cp $SDPBDIR/include/boost/filesystem.hpp  $SDPBDIR/include/boost/filesystem.hpp.old
sed -i '14 a  #include <boost/filesystem/fstream.hpp>' $SDPBDIR/include/boost/filesystem.hpp
```

to fix the headers needed by SPDB IO (`fstream.hpp` was included by default till [266e1ac](https://github.com/boostorg/filesystem/commit/266e1ac892a6f54d807fb35bf639a9aa1c8b2db1)).

### OpenBLAS

Version [0.3.20](https://github.com/xianyi/OpenBLAS/tree/v0.3.20)

```bash
make TARGET=SKYLAKEX
make PREFIX=$SDPBDIR install
```

### LibXML2

Version [2.9.14](https://gitlab.gnome.org/GNOME/libxml2/-/releases/v2.9.14)

```bash
CFLAGS=-O2 ./configure --prefix=$SDPBDIR
make

# optional tests
wget http://www.w3.org/XML/Test/xmlts20080827.tar.gz
tar xvf xmlts20080827.tar.gz
make check

make install
```

### libarchive

Version [3.6.1](https://github.com/libarchive/libarchive/releases/tag/v3.6.1). Standard cmake workflow:

```bash
cmake . -DCMAKE_INSTALL_PREFIX=$SDPBDIR
make

# optional tests
make test

make install
```

### RapidJSON

Version [1.1.0](https://github.com/Tencent/rapidjson/releases/tag/v1.1.0)

```bash
cp -r include $SDPBDIR
chmod 664 $SDPBDIR/include/rapidjson/**.h
chmod 664 $SDPBDIR/include/rapidjson/**/*.h
```

### GMP

Version [6.2.1](https://gmplib.org/download/gmp/gmp-6.2.1.tar.bz2)

```bash
 ./configure --prefix=$SDPBDIR --enable-cxx
 make
 # optional tests
 make check
 make install
```

### MPFR

Version [4.1.0](https://www.mpfr.org/mpfr-4.1.0/). Some tweaks are needed to avoid linking to system-wide libgmp.so (v6.0.0/gcc4.8.5)

```bash
./configure --prefix=$SDPBDIR C_INCLUDE_PATH=$SDPBDIR/include LIBRARY_PATH=$SDPBDIR/lib LD_LIBRARY_PATH=$SDPBDIR/lib --with-gmp=$SDPBDIR
make
# references getting shaky, so tests are very recommended
make check
# if tversion fails it's no bueno
make install
```

### Elemental

Bootstrap collaboration repository, HEAD at [cc51ddb0](https://gitlab.com/bootstrapcollaboration/elemental/-/commit/cc51ddb0f33a467fbf580989d9951620ff2a7877)

```bash
# git checkout cc51ddb0 
mkdir build
cd build
cmake .. -DCMAKE_INSTALL_PREFIX=$SDPBDIR
make
make install
```

## SDPB

Upstream version at the moment of building, HEAD at [4668f72](https://github.com/davidsd/sdpb/commit/4668f72c935e7feba705dd8247d9aacb23185f1c).

Note unusual installation path of libxml2.

```bash
# git checkout 4668f72
./waf configure --prefix=$SDPBDIR \
                --boost-dir=$SDPBDIR \
                --gmpxx-dir=$SDPBDIR \
                --mpfr-dir=$SDPBDIR \
                --elemental-dir=$SDPBDIR \
                --libxml2-dir=$SDPBDIR \
                --rapidjson-dir=$SDPBDIR \
                --libarchive-dir=$SDPBDIR \
                --libxml2-incdir=$SDPBDIR/include/libxml2 \
                --libxml2-libdir=$SDPBDIR/lib
srun ./waf # IO behaves funny when running waf on frontend
./waf install
```

## Usage

MVAPICH2 must be loaded

```bash
module load gcc/8.4.0 mvapich2/2.3.4
```

You can either refer to binaries in the folder, or, for convenience, add them to path.

```bash
export PATH=/work/fsl/sdpb-local/bin:$PATH
```

One may place these two lines in `.bashrc` in the home directory to avoid typing that every time.