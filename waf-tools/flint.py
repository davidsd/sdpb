#! /usr/bin/env python
# encoding: utf-8

def configure(conf):
    # Find Libarchive
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

    if not conf.check_cxx(msg="Checking for flint",
                          fragment='''
#include "flint/flint.h"
#include "flint/arb.h"

int main()
{
    arb_t x;
    arb_init(x);
    arb_const_pi(x, 50 * 3.33);
    arb_printn(x, 50, 0); flint_printf("");
    flint_printf("Computed with FLINT-%s", flint_version);
    arb_clear(x);
}''',
                          includes=flint_incdir,
                          uselib_store='flint',
                          libpath=flint_libdir,
                          rpath=flint_libdir,
                          lib=flint_libs,
                          use=['cxx17']):
        conf.fatal("Could not find flint")


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
