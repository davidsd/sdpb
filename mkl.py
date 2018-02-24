#! /usr/bin/env python
# encoding: utf-8

def configure(conf):
    def get_param(varname,default):
        return getattr(Options.options,varname,'')or default

    mkl_incdir=[]
    mkl_libdir=[]
    mkl_cxxflags=[]
    mkl_linkflags=[]
    mkl_libs=[]

    # Find MKL
    if conf.options.mkl_dir:
        if not conf.options.mkl_incdir:
            mkl_incdir=conf.options.mkl_dir + "/include"
        if not conf.options.mkl_libdir:
            mkl_libdir=conf.options.mkl_dir + "/lib/intel64"

    if conf.options.mkl_incdir:
        mkl_incdir=[conf.options.mkl_incdir]
    if conf.options.mkl_cxxflags:
        mkl_cxxflags=conf.options.mkl_cxxflags.split()

    if conf.options.mkl_libdir:
        mkl_libdir=[conf.options.mkl_libdir]
    if conf.options.mkl_libs:
        mkl_libs=conf.options.mkl_libs.split()
    if not conf.options.mkl_linkflags:
        mkl_linkflags=conf.options.linkflags.split()

    conf.check_cxx(msg="Checking for Mkl",
                   header_name='mkl.h',
                   includes=mkl_incdir,
                   cxxflags=mkl_cxxflags,
                   linkflags=mkl_linkflags,
                   uselib_store='mkl',
                   libpath=mkl_libdir,
                   rpath=mkl_libdir,
                   lib=mkl_libs)

def options(opt):
    mkl=opt.add_option_group('Mkl Options')
    mkl.add_option('--mkl-dir', default='',
                   help='Base directory where mkl is installed')
    mkl.add_option('--mkl-incdir', default='',
                   help='Directory where mkl include files are installed')
    mkl.add_option('--mkl-cxxflags', default='-DMKL_ILP64 -m64',
                   help='Additional flags when compiling (e.g. -DMKL_ILP64 -m64)')
    mkl.add_option('--mkl-libdir', default='',
                   help='Directory where mkl library files are installed')
    mkl.add_option('--mkl-libs', default='mkl_intel_ilp64 mkl_gnu_thread mkl_core gomp pthread m dl',
                   help='Names of the mkl libraries without prefix or suffix (e.g. mkl_intel_ilp64 mkl_gnu_thread mkl_core gomp pthread m dl)')
    mkl.add_option('--mkl-linkflags', default='-Wl,--no-as-needed',
                   help='Additional flags when linking (e.g. -Wl,--no-as-needed)')
    
