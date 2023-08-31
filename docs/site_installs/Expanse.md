# System info

    Rocky Linux release 8.8 (Green Obsidian)
    Linux exp-1-33 4.18.0-477.15.1.el8_8.x86_64 #1 SMP Wed Jun 28 15:04:18 UTC 2023 x86_64 x86_64 x86_64 GNU/Linux

### Useful links

- [Expanse HPC User Guide](https://www.sdsc.edu/support/user_guides/expanse.html)
- [Running MPI jobs](https://www.sdsc.edu/support/user_guides/expanse.html#running)

# Load modules

    module load cpu/0.15.4 gcc/10.2.0 openmpi/4.0.4 gmp/6.1.2 mpfr/4.0.2 cmake/3.18.2 openblas/dynamic/0.3.7 boost/1.74.0 slurm

You may run `module -t list` to view loaded modules,
and `module purge` to unload all modules.

# Install SDPB

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
    ./waf
    ./test/run_all_tests.sh
    ./waf install
    cd ..

# Install scalar_blocks

## Trilinos

    git clone --branch trilinos-release-12-12-branch https://github.com/trilinos/Trilinos.git
    cd Trilinos
    mkdir build
    cd build
    cmake .. -DCMAKE_CXX_COMPILER=g++ -DCMAKE_C_COMPILER=gcc -DCMAKE_Fortran_COMPILER=gfortran -DTrilinos_ENABLE_Sacado=ON -DTrilinos_ENABLE_Kokkos=OFF -DTrilinos_ENABLE_Teuchos=OFF -DCMAKE_INSTALL_PREFIX=$HOME/install
    make && make install
    cd ..

## scalar_blocks

    git clone https://gitlab.com/bootstrapcollaboration/scalar_blocks
    cd scalar_blocks
    ./waf configure --boost-dir=$HOME/install --gmp-dir=$HOME/install --trilinos-dir=$HOME/install --eigen-incdir=/cm/shared/apps/spack/cpu/opt/spack/linux-centos8-zen2/intel-19.1.1.217/eigen-3.3.7-plaog3szjnn3gh6wq5co55xxjuswwo7f/include/eigen3 --prefix=$HOME/install
    ./waf
    ./waf install
    cd ..

# Install blocks_3d

## fmt

    wget https://github.com/fmtlib/fmt/releases/download/6.2.1/fmt-6.2.1.zip
    unzip fmt-6.2.1.zip -d fmt-6.2.1
    cd fmt-6.2.1
    mkdir build
    cd build
    cmake .. -DCMAKE_CXX_COMPILER=g++ -DCMAKE_INSTALL_PREFIX=$HOME/install
    make && make install
    cd ..

## blocks_3d

    git clone https://gitlab.com/bootstrapcollaboration/blocks_3d.git
    cd blocks_3d
    ./waf configure --prefix=$HOME/install --fmt-dir=$HOME/install --fmt-libdir=$HOME/install/lib64 --boost-dir=$HOME/install --eigen-incdir=/cm/shared/apps/spack/cpu/opt/spack/linux-centos8-zen2/intel-19.1.1.217/eigen-3.3.7-plaog3szjnn3gh6wq5co55xxjuswwo7f/include/eigen3
    ./waf
    ./waf install
    cd ..

## Batch Scripts

    /home/wlandry/gnu_openmpi/runs/TTTT_small.sh
    /home/wlandry/gnu_openmpi/runs/blocks.sh
