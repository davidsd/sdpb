#! /usr/bin/env python
# encoding: utf-8

def configure(conf):
    def get_param(varname,default):
        return getattr(Options.options,varname,'')or default

    mkl_incdir=[]
    mkl_libdir=[]
    mkl_cxxflags=['-DMKL_ILP64 -m64']
    mkl_linkflags=['-Wl,--no-as-needed']
    mkl_libs=['mkl_intel_ilp64','mkl_gnu_thread','mkl_core','gomp','pthread','m','dl']

    # Find MKL
    if conf.options.cblas_dir:
        if not conf.options.cblas_incdir:
            mkl_incdir=conf.options.cblas_dir + "/include"
        if not conf.options.cblas_libdir:
            mkl_libdir=conf.options.cblas_dir + "/lib/intel64"

    if conf.options.cblas_incdir:
        mkl_incdir=[conf.options.cblas_incdir]
    if conf.options.cblas_cxxflags:
        mkl_cxxflags=conf.options.cblas_cxxflags.split()

    if conf.options.cblas_libdir:
        mkl_libdir=[conf.options.cblas_libdir]
    if conf.options.cblas_libs:
        mkl_libs=conf.options.cblas_libs.split()
    if conf.options.cblas_linkflags:
        mkl_linkflags=conf.options.cblas_linkflags.split()

    conf.check_cxx(msg="Checking for Mkl",
                   header_name='mkl.h',
                   includes=mkl_incdir,
                   cxxflags=mkl_cxxflags,
                   linkflags=mkl_linkflags,
                   uselib_store='cblas',
                   libpath=mkl_libdir,
                   rpath=mkl_libdir,
                   lib=mkl_libs)

