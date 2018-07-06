#! /usr/bin/env python
# encoding: utf-8

def configure(conf):
    def get_param(varname,default):
        return getattr(Options.options,varname,'')or default

    conf.load('compiler_cxx cxx14')
    
    # Find Elemental
    if conf.options.elemental_dir:
        if not conf.options.elemental_incdir:
            conf.options.elemental_incdir=conf.options.elemental_dir + "/include"
        if not conf.options.elemental_libdir:
            conf.options.elemental_libdir=conf.options.elemental_dir + "/lib"

    if conf.options.elemental_incdir:
        elemental_incdir=conf.options.elemental_incdir.split()
    else:
        elemental_incdir=[]
    if conf.options.elemental_libdir:
        elemental_libdir=conf.options.elemental_libdir.split()
    else:
        elemental_libdir=[]

    if conf.options.elemental_libs:
        elemental_libs=conf.options.elemental_libs.split()
    else:
        elemental_libs=['El', 'pmrrr', 'ElSuiteSparse', 'pthread', 'm', 'mpc',
                        'mpfr', 'gmp', 'metis' ]

    conf.check_cxx(msg="Checking for Elemental",
                   header_name='El.hpp',
                   includes=elemental_incdir,
                   uselib_store='elemental',
                   libpath=elemental_libdir,
                   rpath=elemental_libdir,
                   lib=elemental_libs,
                   use='cxx14')

def options(opt):
    elemental=opt.add_option_group('Elemental Options')
    elemental.add_option('--elemental-dir',
                   help='Base directory where elemental is installed')
    elemental.add_option('--elemental-incdir',
                   help='Directory where elemental include files are installed')
    elemental.add_option('--elemental-libdir',
                   help='Directory where elemental library files are installed')
    elemental.add_option('--elemental-libs',
                   help='Names of the elemental libraries without prefix or suffix\n'
                   '(e.g. "El pmrrr ElSuiteSparse")')
