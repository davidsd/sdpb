#! /usr/bin/env python
# encoding: utf-8

def configure(conf):
    # Find openblas
    import os
    if not conf.options.openblas_incdir:
        for d in ['OPENBLAS_INCLUDE', 'OPENBLAS_INCLUDE_DIR', 'OPENBLAS_INC_DIR']:
            env_dir = os.getenv(d)
            if env_dir:
                conf.to_log('Setting openblas_incdir using environment variable: ' + d + '=' + env_dir)
                conf.options.openblas_incdir = env_dir

    if not conf.options.openblas_libdir:
        for d in ['FLINT_LIB', 'FLINT_LIB_DIR']:
            env_dir = os.getenv(d)
            if env_dir:
                conf.to_log('Setting openblas_libdir using environment variable: ' + d + '=' + env_dir)
                conf.options.openblas_libdir = env_dir

    if conf.options.openblas_dir:
        if not conf.options.openblas_incdir:
            conf.options.openblas_incdir = conf.options.openblas_dir + "/include"
        if not conf.options.openblas_libdir:
            conf.options.openblas_libdir = conf.options.openblas_dir + "/lib"

    if conf.options.openblas_incdir:
        openblas_incdir = conf.options.openblas_incdir.split()
    else:
        openblas_incdir = []
    if conf.options.openblas_libdir:
        openblas_libdir = conf.options.openblas_libdir.split()
    else:
        openblas_libdir = []

    if conf.options.openblas_libs:
        openblas_libs = conf.options.openblas_libs.split()
    else:
        openblas_libs = ['openblas']

    if not conf.check_cxx(msg="Checking for openblas",
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
                          includes=openblas_incdir,
                          uselib_store='openblas',
                          libpath=openblas_libdir,
                          rpath=openblas_libdir,
                          lib=openblas_libs,
                          use=['cxx17']):
        conf.fatal("Could not find openblas")


def options(opt):
    openblas = opt.add_option_group('openblas Options')
    openblas.add_option('--openblas-dir',
                        help='Base directory where openblas is installed')
    openblas.add_option('--openblas-incdir',
                        help='Directory where openblas include files are installed')
    openblas.add_option('--openblas-libdir',
                        help='Directory where openblas library files are installed')
    openblas.add_option('--openblas-libs',
                        help='Names of the openblas and openblas libraries without prefix or suffix\n'
                             '(e.g. "openblas")')
