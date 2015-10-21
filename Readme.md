## Contents

* [SDPB](#sdpb)
* [Installation and Requirements](#installation-and-requirements)
  * [Linux Installation](#linux-installation)
  * [Mac OS X Installation](#mac-os-x-installation)
  * [Windows Installation](#windows-installation)
* [Usage](#usage)
* [Changelog](#changelog)
* [Attribution](#attribution)
* [Acknowledgements](#acknowledgements)
* [Works Using SDPB](#works-using-sdpb)

# SDPB

SDPB is an open-source, arbitrary-precision, parallelized semidefinite
program solver, designed for the conformal bootstrap.

It solves the following type of optimization problem

```
maximize:  b_0 + \sum_n b_n y_n over (y_1,...,y_N),

such that: M_{0j}(x) + \sum_n y_n M_{nj}(x) is positive semidefinite
           for all x >= 0 and 1 <= j <= J,

where each M_{nj}(x) is a polynomial matrix in x.
```

For more information, see [A Semidefinite Program Solver for the Conformal Bootstrap](http://arxiv.org/abs/1502.02033)
and [the manual](https://github.com/davidsd/sdpb/blob/master/docs/SDPB-Manual.pdf).

Author: David Simmons-Duffin (davidsd@gmail.com). As of February 2015, I am
supported by DOE grant number DE-SC0009988 and a William D. Loughlin Membership
at the [Institute for Advanced Study](http://sns.ias.edu).

## Installation and Requirements

SDPB requires

- [Boost C++ Libraries](http://www.boost.org/) (version 1.54 or later).

- [The GNU Multiprecision Library](https://gmplib.org/).

System-specific installation instructions are below.  Please use the
issue tracker if you encounter installation problems. For those with
experience packaging software, I'd appreciate help making SDPB easier
to install.

### Linux Installation

SDPB has been tested on Red Hat Linux. (Many thanks to Chris Beem for helping with installation instructions.) To install,

1. Download Boost and GMP from the links above. Install GMP with the option `--enable-cxx` added to `./configure`. Install Boost.

2. Edit the `Makefile` to define the variables `GMPINCLUDEDIR`,
`BOOSTINCLUDEDIR`, and `LIBDIR.` Ensure `LIBDIR` is in your `LD_LIBRARY_PATH`. 

3. Type `make` to build the `sdpb` executable.

### Mac OS X Installation

The following instructions have been tested on Mac OS 10.10 Yosemite.  (Many thanks to Ying Lin.)

1. Install Homebrew and `gcc-4.9` (or later), for instance by running the following commands

        # Install Homebrew
        ruby -e "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/master/install)"

        # Update Ruby
        brew install ruby

        # Install the latest gcc as well as its dependencies
        # the option --without-multilib avoids a bug in OpenMP support
        brew install gcc --without-multilib

2. Make `/usr/local/bin/g++-4.9` (or whatever version you have) the default compiler by renaming `gcc` and `g++` in `/usr/bin` and creating symlinks

        ln -s /usr/local/bin/g++-4.9 /usr/local/bin/g++
        ln -s /usr/local/bin/gcc-4.9 /usr/local/bin/gcc

3. Unfortunately, Homebrew's versions of GMP and Boost will not work -- they must be compiled from source. Download the latest GMP from [the GMP website](https://gmplib.org/). Upack the tarball (you may need `lzip` which you can install with `brew install lzip`) and `cd` to the `gmp` directory.  Run

        ./configure --enable-cxx
        make
        make install

4. Download Boost from [the Boost website](http://www.boost.org/).  Unpack the tarball and `cd` to the `boost` directory. Run

        ./bootstrap.sh
        ./b2
        sudo ./b2 install
        
   (Note that `bootstrap.sh` above is just an installation script and has absolutely nothing
   to do with the conformal bootstrap -- lots of people like the name "bootstrap"!)
        
5. Type `make` in the `sdpb` directory to compile the `sdpb` executable. If you installed any of the
above software in custom locations, you'll have to modify variables in the
`Makefile` as described in the Linux instructions.

### Windows Installation

The following instructions have been tested on Windows 8.1 using Boost 1.59.0, GMP 6.0.0, Cygwin DLL ver. 2.2.1. These are written to be accessible to people with no Unix experience by a person with little Unix experience, so suggestions to improve are appreciated.

Below are the steps to build SDPB executable on your machine. These steps use Cygwin to provide POSIX environment for SDPB and its dependecies. To simplify the process, it is assumed that you do not intend to use Boost or GMP installations for a purpose other than building SDPB. The goal here was not to follow the best practices of Unix, but to make the installation process as quick as possible.

It is not required to have any specific directory structure for the installation. However, to be concrete, it is assumed in the instructions that you have created a directory `C:\SDPB`, where you are going to store files related to SDPB. The filenames used in the instructions are suitable for Boost 1.59.0 and GMP 6.0.0a. If you are using newer versions, you will have to use the appropriate file and directory names.

1. Download and install Cygwin from [the Cygwin installation webpage](http://cygwin.com/install.html). Choose the appropriate version (32-bit or 64-bit) and run the dowloaded executable file (`setup-x86.exe` or `setup-x86_64.exe`). 
  1. Choose `Install from Internet`.
  2. Choose the installation directory and the directory for storing the downloaded packages.
  3. Choose appropriate proxy settings. If you don't know what it is, try `Direct Connection` or `Use Internet Explorer Proxy Settings`.
  4. Choose a download site. Any choice should work.
  5. Choose the packages for installation. You will need
    * All the `Base` packages (choosen by defaul, no action should be required).
    * `bzip2` and `lzip` in `Archive`.
    * `gcc-core`, `gcc-g++` and `make` under `Devel`
    * `m4` under `Interpreters`
    * If in further steps some required tools are found to be missing, they can be installed by re-running this installer.
  6. Continue and accept any package dependencies.
  7. (optional) When the installation is finished, add the `\bin` subdirectory of the Cygwin installation directory to your system `Path` variable.
2. Download and build Boost from [the Boost website](http://www.boost.org). You should get the Unix variant. Put the downloaded file in `C:\SDPB`.
  1. Run the Cygwin terminal. Navigate to `C:\SDPB` by typing

           cd /cygdrive/c/SDPB

     and pressing Enter key. Note that the paths are case-sensitive, but the drive letter has to be lowercase.
  2. Unpack the Boost archive by typing

           tar --bzip2 -xf boost_1_59_0.tar.bz2

  3. Navigate to Boost directory and build the required Boost libraries. For this, type

           cd boost_1_59_0/
           ./bootstrap.sh --with-libraries=filesystem,serialization,program_options,date_time,timer
           ./b2 stage
           
     Note: Using `stage` target instead of `install` can save a lot of time by skipping the copying of Boost header files, which is done more efficiently by Windows methods.
           
3. Dowload and build GMP from [the GMP website](https://gmplib.org). Download the latest version, and put the file in `C:\SDPB`.
  1. In Cygwin terminal, navigate to `C:\SDPB` and upack the archive. For this, type

           cd /cygdrive/c/SDPB
           tar --lzip -xf gmp-6.0.0a.tar.lz
           
  2. Build GMP libraries. For this, type

           mkdir installation
           cd gmp-6.0.0
           ./configure --enable-cxx --prefix=/cygdrive/c/SDPB/installation
           make install
           
  3. Run GMP tests to make sure that everything is in order

           make check
           
4. Collect the header and lib files in one place
  1. Move `C:\SDPB\boost_1_59_0\boost` directory into `C:\SDPB\installation\include` (so you will now have `C:\SDPB\installation\include\boost` directory).
  2. Copy the contents of `C:\SDPB\boost_1_59_0\stage\lib` directory into `C:\SDPB\installation\lib`.
  3. (optional) Add `C:\SDPB\installation\lib` to your system `Path` variable.

5. Download and build SDPB. Obtain the latest SDPB sources from github, for example, by dowloading the zip of the repository and unpacking the contents into `C:\SDPB`. As the result you should have the contents of the repository in `C:\SDPB\sdpb-master`.
  1. Use any text-editing software to open `C:\SDPB\sdpb-master\Makefile`
    * Follow the instructions in the file to locate the variables `GMPINCLUDEDIR` and `BOOSTINCLUDEDIR` as well as `LIBDIR`, which should be edited on a **Windows** system.
    * Edit the appropriate lines to look as follows

                      GMPINCLUDEDIR = /cygdrive/c/SDPB/installation/include
                      BOOSTINCLUDEDIR = /cygdrive/c/SDPB/installation/include/boost
                      LIBDIR = /cygdrive/c/SDPB/installation/lib

    * Save the file and exit
  2. In Cygwin terminal, navigate to `C:\SDPB\sdpb-master` and build SDPB
  
           cd /cygdrive/c/SDPB/sdpb-master
           make

    After the process finishes, `sdpb.exe` should appear in `C:\SDPB\sdpb-master`

6. If you have added all the suggested locations to your system `Path` variable, then you should be able to run sdpb.exe. If you have not, then you need to either modify `Path` variable or put copies of the `.dll` files from `C:\SDPB\installation\lib` and a copy of `cygwin1.dll` from `\bin` subdirectory of Cygwin installation in the directory where you would like to keep `sdpb.exe`. You then might want to add this directory to `Path` variable.



## Usage

Type `sdpb --help` for the syntax and a list of options.

The input format for SDPB is XML-based and described in the manual.
The Mathematica file `SDPB.m` includes code to export semidefinite
programs in this format, along with some examples.

## Changelog

- 5/29/15: Fixed a bug in `SDP.cpp` that caused SDPB to incorrectly load
  matrices with dimensions larger than 2x2. (Thanks to Filip Kos for
  the fix.)
- 10/20/15: Fixed a Windows incompatibility in the parser ([Pull request](https://github.com/davidsd/sdpb/pull/8))

## Attribution

If you use SDPB in work that results in publication, please cite

- D. Simmons-Duffin, "A Semidefinite Program Solver for the
  Conformal Bootstrap", [arXiv:1502.02033 \[hep-th\]](http://arxiv.org/abs/1502.02033).

Depending on how SDPB is used, please also consider the following sources:

The first use of semidefinite programming in the bootstrap:

- D. Poland, D. Simmons-Duffin and A. Vichi, "Carving Out the Space of
  4D CFTs," JHEP 1205, 110 (2012) [arXiv:1109.5176 \[hep-th\]](http://arxiv.org/abs/1109.5176).

The generalization of semidefinite programming methods to arbitrary
spacetime dimension:

- F. Kos, D. Poland and D. Simmons-Duffin, "Bootstrapping the O(N)
  Vector Models," JHEP 1406, 091 (2014) [arXiv:1307.6856 \[hep-th\]](http://arxiv.org/abs/1307.6856).

The generalization of semidefinite programming methods to arbitrary
systems of correlation functions:

- F. Kos, D. Poland and D. Simmons-Duffin, "Bootstrapping Mixed
  Correlators in the 3D Ising Model," [arXiv:1406.4858 \[hep-th\]](http://arxiv.org/abs/1406.4858).

## Acknowledgements

- SDPB Makes extensive use of [MPACK](http://mplapack.sourceforge.net/), the multiple precision linear algebra library written by Nakata Maho.  Several source files from MPACK are included in the SDPB source tree (see the license at the top of those files).

- SDPB uses Lee Thomason's [tinyxml2](http://www.grinninglizard.com/tinyxml2/) library (included in the source tree) for parsing.

- The design of SDPB was partially based on the solvers [SDPA](http://sdpa.sourceforge.net/) and SDPA-GMP, which were essential sources of inspiration and examples.

- Thanks to Filip Kos, David Poland, and Alessandro Vichi for collaboration in developing semidefinite programming methods for the conformal bootstrap and assistance testing SDPB.

- Thanks to Amir Ali Ahmadi, Hande Benson, Pablo Parrilo, and Robert Vanderbei for advice and discussions about semidefinite programming.

- Thanks also to Noah Stein, who first suggested the idea of semidefinite programming to me in [this Math Overflow question](http://mathoverflow.net/questions/33242/continuous-linear-programming-estimating-a-solution).

## Works Using SDPB

- F. Kos, D. Poland, D. Simmons-Duffin and A. Vichi,
  "Bootstrapping the O(N) Archipelago,"
  [arXiv:1504.07997 [hep-th]](http://arxiv.org/abs/1504.07997).

- S. M. Chester, S. Giombi, L. V. Iliesiu, I. R. Klebanov, S. S. Pufu and R. Yacoby,
  "Accidental Symmetries and the Conformal Bootstrap,"
  [arXiv:1507.04424 [hep-th]](http://arxiv.org/abs/1507.04424).

- C. Beem, M. Lemos, L. Rastelli and B. C. van Rees,
  "The (2,0) superconformal bootstrap,"
  [arXiv:1507.05637 [hep-th]](http://arxiv.org/abs/1507.05637).

- L. Iliesiu, F. Kos, D. Poland, S. S. Pufu, D. Simmons-Duffin and R. Yacoby,
  "Bootstrapping 3D Fermions,"
  [arXiv:1508.00012 [hep-th]](http://arxiv.org/abs/1508.00012).

- D. Poland and A. Stergiou,
  "Exploring the Minimal 4D N=1 SCFT,"
  [arXiv:1509.06368 [hep-th]](http://arxiv.org/abs/1509.06368).

-  M. Lemos and P. Liendo,
  "Bootstrapping N=2 chiral correlators,"
  [arXiv:1510.03866 [hep-th]](http://arxiv.org/abs/1510.03866).

(Please let me know if you would like your paper to be included in this list!)
