#! /usr/bin/env python
# encoding: utf-8

from check_config import check_config


def configure(conf):
    # Find FLINT
    import os

    if not conf.options.flint_incdir:
        for d in ['FLINT_INCLUDE', 'FLINT_INCLUDE_DIR', 'FLINT_INC_DIR']:
            env_dir = os.getenv(d)
            if env_dir:
                conf.to_log('Setting flint_incdir using environment variable: ' + d + '=' + env_dir)
                conf.options.flint_incdir = env_dir

    if not conf.options.flint_libdir:
        for d in ['FLINT_LIB', 'FLINT_LIB_DIR']:
            env_dir = os.getenv(d)
            if env_dir:
                conf.to_log('Setting flint_libdir using environment variable: ' + d + '=' + env_dir)
                conf.options.flint_libdir = env_dir

    if conf.options.flint_dir:
        if not conf.options.flint_incdir:
            conf.options.flint_incdir = conf.options.flint_dir + "/include"
        if not conf.options.flint_libdir:
            conf.options.flint_libdir = conf.options.flint_dir + "/lib"
            conf.options.flint_libdir += " " + conf.options.flint_dir + "/lib64"

    if conf.options.flint_incdir:
        flint_incdir = conf.options.flint_incdir.split()
    else:
        flint_incdir = []
    if conf.options.flint_libdir:
        flint_libdir = conf.options.flint_libdir.split()
    else:
        flint_libdir = []

    if conf.options.flint_libs:
        flint_libs = conf.options.flint_libs.split()
    else:
        flint_libs = ['flint']

    check_config(conf,
                 fragment='''
#include <gmp.h>
#include "flint/flint.h"
#include "flint/fmpz.h"
#include <vector>

int main()
{
    static_assert(__FLINT_RELEASE >= 20800, "FLINT 2.8.0 or later required, current version: " FLINT_VERSION);

    std::vector<mp_limb_t> primes{2, 3, 5};
    fmpz_comb_t comb;
    fmpz_comb_init(comb, primes.data(), primes.size());
    flint_printf("Computed with FLINT-%s", FLINT_VERSION);
}''',
                 includes=flint_incdir,
                 uselib_store='flint',
                 libpath=flint_libdir,
                 lib=flint_libs,
                 use=['cxx17', 'gmpxx', 'mpfr', 'cblas'])


def options(opt):
    flint = opt.add_option_group('FLINT Options')
    flint.add_option('--flint-dir',
                     help='Base directory where flint is installed')
    flint.add_option('--flint-incdir',
                     help='Directory where flint include files are installed')
    flint.add_option('--flint-libdir',
                     help='Directory where flint library files are installed')
    flint.add_option('--flint-libs',
                     help='Names of the flint and flint libraries without prefix or suffix\n'
                          '(e.g. "flint")')
