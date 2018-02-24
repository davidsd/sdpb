#! /usr/bin/env python
# encoding: utf-8

def configure(conf):
    def get_param(varname,default):
        return getattr(Options.options,varname,'')or default

    # Find CBLAS
    if conf.options.cblas_dir:
        if not conf.options.cblas_incdir:
            conf.options.cblas_incdir=conf.options.cblas_dir + "/include"
        if not conf.options.cblas_libdir:
            conf.options.cblas_libdir=conf.options.cblas_dir + "/lib"

    if conf.options.cblas_incdir:
        cblas_incdir=[conf.options.cblas_incdir]
    else:
        cblas_incdir=[]
    if conf.options.cblas_libdir:
        cblas_libdir=[conf.options.cblas_libdir]
    else:
        cblas_libdir=[]

    if conf.options.cblas_libs:
        cblas_libs=conf.options.cblas_libs.split()
    else:
        cblas_libs=['blas']

    conf.check_cxx(msg="Checking for Cblas",
                   header_name='cblas.h',
                   includes=cblas_incdir,
                   uselib_store='cblas',
                   libpath=cblas_libdir,
                   rpath=cblas_libdir,
                   lib=cblas_libs)

def options(opt):
    cblas=opt.add_option_group('Cblas Options')
    cblas.add_option('--cblas-dir',
                   help='Base directory where cblas is installed')
    cblas.add_option('--cblas-incdir',
                   help='Directory where cblas include files are installed')
    cblas.add_option('--cblas-libdir',
                   help='Directory where cblas library files are installed')
    cblas.add_option('--cblas-libs',
                   help='Names of the cblas libraries without prefix or suffix\n'
                   '(e.g. "cblas")')
