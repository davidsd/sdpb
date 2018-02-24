#! /usr/bin/env python
# encoding: utf-8

def configure(conf):
    def get_param(varname,default):
        return getattr(Options.options,varname,'')or default

    # Find MKL
    if conf.options.mkl_dir:
        if not conf.options.mkl_incdir:
            conf.options.mkl_incdir=conf.options.mkl_dir + "/include"
        if not conf.options.mkl_libdir:
            conf.options.mkl_libdir=conf.options.mkl_dir + "/lib/intel64"

    if conf.options.mkl_incdir:
        mkl_incdir=[conf.options.mkl_incdir]
    else:
        mkl_incdir=[]
    if conf.options.mkl_libdir:
        mkl_libdir=[conf.options.mkl_libdir]
    else:
        mkl_libdir=[]

    if conf.options.mkl_libs:
        mkl_libs=conf.options.mkl_libs.split()
    else:
        mkl_libs=['mkl_intel_ilp64', 'mkl_gnu_thread', 'mkl_core', 'gomp', 'pthread',
                  'm', 'dl']

    conf.check_cxx(msg="Checking for Mkl",
                   header_name='mkl.h',
                   includes=mkl_incdir,
                   cxxflags=['-DMKL_ILP64', '-m64'],
                   linkflags=['-Wl,--no-as-needed'],
                   uselib_store='mkl',
                   libpath=mkl_libdir,
                   rpath=mkl_libdir,
                   lib=mkl_libs)

def options(opt):
    mkl=opt.add_option_group('Mkl Options')
    mkl.add_option('--mkl-dir',
                   help='Base directory where mkl is installed')
    mkl.add_option('--mkl-incdir',
                   help='Directory where mkl include files are installed')
    mkl.add_option('--mkl-libdir',
                   help='Directory where mkl library files are installed')
    mkl.add_option('--mkl-libs',
                   help='Names of the mkl libraries without prefix or suffix\n'
                   '(e.g. "mkl_intel_ilp64 mkl_gnu_thread mkl_core gomp pthread m dl")')

#     Link: -L${MKLROOT}/lib/intel64 -Wl,--no-as-needed -lmkl_intel_ilp64 -lmkl_gnu_thread -lmkl_core -lgomp -lpthread -lm -ldl

# Compile:  -DMKL_ILP64 -m64 -I${MKLROOT}/include
