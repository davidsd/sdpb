#! /usr/bin/env python
# encoding: utf-8

from check_config import check_config


def configure(conf):
    import os

    # Find Elemental
    if conf.options.elemental_dir:
        if not conf.options.elemental_incdir:
            conf.options.elemental_incdir = conf.options.elemental_dir + "/include"
        if not conf.options.elemental_libdir:
            lib = conf.options.elemental_dir + "/lib"
            if os.path.isdir(lib):
                conf.options.elemental_libdir = lib
            lib64 = conf.options.elemental_dir + "/lib64"
            if os.path.isdir(lib64):
                if conf.options.elemental_libdir:
                    conf.options.elemental_libdir += " " + lib64
                else:
                    conf.options.elemental_libdir = lib64

    if conf.options.elemental_incdir:
        elemental_incdir = conf.options.elemental_incdir.split()
    else:
        elemental_incdir = []
    if conf.options.elemental_libdir:
        elemental_libdir = conf.options.elemental_libdir.split()
    else:
        elemental_libdir = []

    if conf.options.elemental_libs:
        elemental_libs = conf.options.elemental_libs.split()
    else:
        elemental_libs = ['El', 'pmrrr', 'ElSuiteSparse']

    check_config(conf,
                 fragment="#include <El.hpp>\nint main(int argc, char* argv[]) {El::Environment env( argc, argv ); El::BigFloat big;}\n",
                 includes=elemental_incdir,
                 uselib_store='elemental',
                 libpath=elemental_libdir,
                 lib=elemental_libs,
                 use=['cxx17', 'gmpxx', 'mpfr'])


def options(opt):
    elemental = opt.add_option_group('Elemental Options')
    elemental.add_option('--elemental-dir',
                         help='Base directory where elemental is installed')
    elemental.add_option('--elemental-incdir',
                         help='Directory where elemental include files are installed')
    elemental.add_option('--elemental-libdir',
                         help='Directory where elemental library files are installed')
    elemental.add_option('--elemental-libs',
                         help='Names of the elemental libraries without prefix or suffix\n'
                              '(e.g. "El pmrrr ElSuiteSparse")')
