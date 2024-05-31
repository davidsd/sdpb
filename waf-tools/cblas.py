#! /usr/bin/env python
# encoding: utf-8

from check_config import check_config


def configure(conf):
    # Find cblas
    import os

    if not conf.options.cblas_incdir:
        for libname in ['CBLAS', 'OPENBLAS']:
            for suffix in ['_INCLUDE', '_INCLUDE_DIR', '_INC_DIR', '_INCDIR']:
                d = libname + suffix
                env_dir = os.getenv(d)
                if env_dir:
                    conf.to_log('Setting cblas_incdir using environment variable: ' + d + '=' + env_dir)
                    conf.options.cblas_incdir = env_dir

    if not conf.options.cblas_libdir:
        for libname in ['CBLAS', 'OPENBLAS']:
            for suffix in ['_LIB', '_LIB_DIR', '_LIBDIR']:
                d = libname + suffix
                env_dir = os.getenv(d)
                if env_dir:
                    conf.to_log('Setting cblas_libdir using environment variable: ' + d + '=' + env_dir)
                    conf.options.cblas_libdir = env_dir

    if conf.options.cblas_dir:
        if not conf.options.cblas_incdir:
            conf.options.cblas_incdir = conf.options.cblas_dir + "/include"
        if not conf.options.cblas_libdir:
            conf.options.cblas_libdir = ' '.join([conf.options.cblas_dir + "/lib", conf.options.cblas_dir + "/lib64"])

    if conf.options.cblas_incdir:
        cblas_incdir = conf.options.cblas_incdir.split()
    else:
        cblas_incdir = []
    if conf.options.cblas_libdir:
        cblas_libdir = conf.options.cblas_libdir.split()
    else:
        cblas_libdir = []

    if conf.options.cblas_libs:
        cblas_libs = conf.options.cblas_libs.split()
    else:
        cblas_libs = ['openblas']

    check_config(conf,
                 fragment='''
#include "cblas.h"

int main()
{
    int m=10,n=10,k=10;
    double x[m*k];
    double y[k*n];
    double result[m*n];
    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, m, n, k, 1.0, x, k, y, n, 0.0, result, n);
}''',
                 includes=cblas_incdir,
                 uselib_store='cblas',
                 libpath=cblas_libdir,
                 lib=cblas_libs,
                 packages=cblas_libs,
                 use=['cxx17'])


def options(opt):
    cblas = opt.add_option_group('CBLAS Options')
    cblas.add_option('--cblas-dir',
                     help='Base directory where CBLAS implementation (e.g. OpenBLAS) is installed')
    cblas.add_option('--cblas-incdir',
                     help='Directory where CBLAS include files are installed')
    cblas.add_option('--cblas-libdir',
                     help='Directory where CBLAS library files are installed')
    cblas.add_option('--cblas-libs',
                     help='Names of the CBLAS libraries without prefix or suffix\n'
                          '(e.g. "openblas")')
