
* [Requirements](#requirements)
* [Installation](#installation)
  * [Linux](#linux)
  * [Mac OS X](#mac-os-x)
  * [Windows](#windows)

# Requirements

SDPB requires

- [Boost C++ Libraries](http://www.boost.org/) (version 1.54 or later).

- [The GNU Multiprecision Library](https://gmplib.org/).

System-specific installation instructions are below.  Please use the
issue tracker if you encounter installation problems. For those with
experience packaging software, I'd appreciate help making SDPB easier
to install.

# Installation

## Linux

(Thanks to Chris Beem for help with these instructions.)

SDPB has been tested on Red Hat Linux. To install,

1. Download Boost and GMP from the links above. Install GMP with the option `--enable-cxx` added to `./configure`. Install Boost.

2. Edit the `Makefile` to define the variables `GMPINCLUDEDIR`,
`BOOSTINCLUDEDIR`, and `LIBDIR.` Ensure `LIBDIR` is in your `LD_LIBRARY_PATH`. 

3. Type `make` to build the `sdpb` executable.

## Mac OS X

(Thanks to Ying Lin for these instructions.)

The following instructions have been tested on Mac OS 10.10 Yosemite.

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

## Windows

(Thanks to Peter Kravchuk for these instructions.)

The following instructions have been tested on Windows 8.1 using Boost 1.59.0, GMP 6.0.0, Cygwin DLL ver. 2.2.1. These are written to be accessible to people with no Unix experience by a person with little Unix experience, so suggestions to improve are appreciated.

Below are the steps to build SDPB executable on your machine. These steps use Cygwin to provide POSIX environment for SDPB and its dependecies. To simplify the process, it is assumed that you do not intend to use Boost or GMP installations for a purpose other than building SDPB. The goal here was not to follow the best practices of Unix, but to make the installation process as quick as possible.

It is not required to have any specific directory structure for the installation. However, to be concrete, it is assumed in the instructions that you have created a directory `C:\SDPB`, where you are going to store files related to SDPB. The filenames used in the instructions are suitable for Boost 1.59.0 and GMP 6.0.0a. If you are using newer versions, you will have to use the appropriate file and directory names.

1. Download and install Cygwin from [the Cygwin installation webpage](http://cygwin.com/install.html). Choose the appropriate version (32-bit or 64-bit) and run the dowloaded executable file (`setup-x86.exe` or `setup-x86_64.exe`). 
  1. Choose `Install from Internet`.
  2. Choose the installation directory and the directory for storing the downloaded packages.
  3. Choose appropriate proxy settings. If you don't know what it is, try `Direct Connection` or `Use Internet Explorer Proxy Settings`.
  4. Choose a download site. Any choice should work.
  5. Choose the packages for installation. You will need
    * All the `Base` packages (choosen by default, no action should be required).
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
  1. In Cygwin terminal, navigate to `C:\SDPB` and unpack the archive. For this, type

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

5. Download and build SDPB. Obtain the latest SDPB sources from github, for example, by downloading the zip of the repository and unpacking the contents into `C:\SDPB`. As the result you should have the contents of the repository in `C:\SDPB\sdpb-master`.
  1. Use any text-editing software to open `C:\SDPB\sdpb-master\Makefile`
    * Follow the instructions in the file to locate those instances of the variables `GMPINCLUDEDIR`, `BOOSTINCLUDEDIR` and `LIBDIR` which should be edited on a **Windows** system.
    * Edit the appropriate lines to look as follows

                      GMPINCLUDEDIR = /cygdrive/c/SDPB/installation/include
                      BOOSTINCLUDEDIR = /cygdrive/c/SDPB/installation/include/boost
                      LIBDIR = /cygdrive/c/SDPB/installation/lib

    * Save the file and exit
  2. In Cygwin terminal, navigate to `C:\SDPB\sdpb-master` and build SDPB
  
           cd /cygdrive/c/SDPB/sdpb-master
           make

    After the process finishes, `sdpb.exe` should appear in `C:\SDPB\sdpb-master`

6. If you have added all the suggested locations to your system `Path` variable, then you should be able to run `sdpb.exe`. If you have not, then you need to either modify `Path` variable or put copies of the `.dll` files from `C:\SDPB\installation\lib` and a copy of `cygwin1.dll` from `\bin` subdirectory of Cygwin installation into the directory where you would like to keep `sdpb.exe`. You then might want to add this directory to `Path` variable to be able to access SDPB from Windows Command Line at any time.


