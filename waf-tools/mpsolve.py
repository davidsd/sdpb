#! /usr/bin/env python
# encoding: utf-8

from check_config import check_config


def configure(conf):
    import os

    if not conf.options.mpsolve_incdir:
        for d in ['MPSOLVE_INCLUDE', 'MPSOLVE_INCLUDE_DIR', 'MPSOLVE_INC_DIR']:
            env_dir = os.getenv(d)
            if env_dir:
                conf.to_log('Setting mpsolve_incdir using environment variable: ' + d + '=' + env_dir)
                conf.options.mpsolve_incdir = env_dir

    if not conf.options.mpsolve_libdir:
        for d in ['MPSOLVE_LIB', 'MPSOLVE_LIB_DIR']:
            env_dir = os.getenv(d)
            if env_dir:
                conf.to_log('Setting mpsolve_libdir using environment variable: ' + d + '=' + env_dir)
                conf.options.mpsolve_libdir = env_dir

    if not conf.options.mpsolve_dir:
        env_dir = os.getenv('MPSOLVE_DIR')
        if env_dir:
            conf.to_log('Setting mpsolve_dir using environment variable: MPSOLVE_DIR=' + env_dir)
            conf.options.mpsolve_dir = env_dir

    # Find MPSOLVE
    if conf.options.mpsolve_dir:
        conf.to_log('Using mpsolve_dir: ' + conf.options.mpsolve_dir)
        if not conf.options.mpsolve_incdir:
            conf.to_log('Setting mpsolve_incdir using mpsolve_dir: ' + conf.options.mpsolve_dir)
            conf.options.mpsolve_incdir = conf.options.mpsolve_dir + "/include"
        if not conf.options.mpsolve_libdir:
            conf.to_log('Setting mpsolve_libdir using mpsolve_dir: ' + conf.options.mpsolve_dir)
            conf.options.mpsolve_libdir = conf.options.mpsolve_dir + "/lib"

    if conf.options.mpsolve_incdir:
        conf.to_log('Using mpsolve_incdir: ' + conf.options.mpsolve_incdir)
        mpsolve_incdir = conf.options.mpsolve_incdir.split()
    else:
        mpsolve_incdir = []
    if conf.options.mpsolve_libdir:
        conf.to_log('Using mpsolve_libdir: ' + conf.options.mpsolve_libdir)
        mpsolve_libdir = conf.options.mpsolve_libdir.split()
    else:
        mpsolve_libdir = []

    if conf.options.mpsolve_libs:
        conf.to_log('Using mpsolve_libs: ' + conf.options.mpsolve_libs)
        mpsolve_libs = conf.options.mpsolve_libs.split()
    else:
        mpsolve_libs = ['mps']

    check_config(conf,
                 fragment='''
#include <mps/mps.h>
#include <gmp.h>

int main()
{
  long int n = 5;
  mps_context *s = mps_context_new();
  mps_context_select_algorithm(s, MPS_ALGORITHM_STANDARD_MPSOLVE);
  
  // The following functions are missing in MPSolve 3.2.1 and earlier.
  // They have been added to 3.2.2 by the PR:
  // https://github.com/robol/MPSolve/pull/42
  mps_context_set_search_set(s, MPS_SEARCH_SET_POSITIVE_REAL_PART);
  mps_context_set_n_threads(s, 1);
  
  mps_monomial_poly * p = mps_monomial_poly_new(s, n);
  mpq_t one;
  mpq_init(one);
}
                 ''',
                 includes=mpsolve_incdir,
                 uselib_store='mpsolve',
                 libpath=mpsolve_libdir,
                 lib=mpsolve_libs,
                 use=['cxx17', 'gmpxx'])


def options(opt):
    mpsolve = opt.add_option_group('MPSolve Options')
    mpsolve.add_option('--mpsolve-dir',
                       help='Base directory where mpsolve is installed')
    mpsolve.add_option('--mpsolve-incdir',
                       help='Directory where mpsolve include files are installed')
    mpsolve.add_option('--mpsolve-libdir',
                       help='Directory where mpsolve library files are installed')
    mpsolve.add_option('--mpsolve-libs',
                       help='Names of the mpsolve libraries without prefix or suffix\n'
                            '(e.g. "mps")')
