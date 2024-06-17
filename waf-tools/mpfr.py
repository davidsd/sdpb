#! /usr/bin/env python
# encoding: utf-8

from check_config import check_config


def configure(conf):
    import os

    if not conf.options.mpfr_incdir:
        for d in ['MPFR_INCLUDE', 'MPFR_INCLUDE_DIR', 'MPFR_INC_DIR']:
            env_dir = os.getenv(d)
            if env_dir:
                conf.to_log('Setting mpfr_incdir using environment variable: ' + d + '=' + env_dir)
                conf.options.mpfr_incdir = env_dir

    if not conf.options.mpfr_libdir:
        for d in ['MPFR_LIB', 'MPFR_LIB_DIR']:
            env_dir = os.getenv(d)
            if env_dir:
                conf.to_log('Setting mpfr_libdir using environment variable: ' + d + '=' + env_dir)
                conf.options.mpfr_libdir = env_dir

    if not conf.options.mpfr_dir:
        env_dir = os.getenv('MPFR_DIR')
        if env_dir:
            conf.to_log('Setting mpfr_dir using environment variable: ' + d + '=' + env_dir)
            conf.options.mpfr_dir = env_dir

    # Find MPFR
    if conf.options.mpfr_dir:
        conf.to_log('Using mpfr_dir: ' + conf.options.mpfr_dir)
        if not conf.options.mpfr_incdir:
            conf.to_log('Setting mpfr_incdir using mpfr_dir: ' + conf.options.mpfr_dir)
            conf.options.mpfr_incdir = conf.options.mpfr_dir + "/include"
        if not conf.options.mpfr_libdir:
            conf.to_log('Setting mpfr_libdir using mpfr_dir: ' + conf.options.mpfr_dir)
            conf.options.mpfr_libdir = conf.options.mpfr_dir + "/lib"

    if conf.options.mpfr_incdir:
        conf.to_log('Using mpfr_incdir: ' + conf.options.mpfr_incdir)
        mpfr_incdir = conf.options.mpfr_incdir.split()
    else:
        mpfr_incdir = []
    if conf.options.mpfr_libdir:
        conf.to_log('Using mpfr_libdir: ' + conf.options.mpfr_libdir)
        mpfr_libdir = conf.options.mpfr_libdir.split()
    else:
        mpfr_libdir = []

    if conf.options.mpfr_libs:
        conf.to_log('Using mpfr_libs: ' + conf.options.mpfr_libs)
        mpfr_libs = conf.options.mpfr_libs.split()
    else:
        mpfr_libs = ['mpfr']

    check_config(conf,
                 fragment='''
#include <mpfr.h>
int main()
{
  mpfr_t x;
  mpfr_init(x);
}
                 ''',
                 includes=mpfr_incdir,
                 uselib_store='mpfr',
                 libpath=mpfr_libdir,
                 lib=mpfr_libs,
                 use=['cxx17', 'gmpxx'])


def options(opt):
    mpfr = opt.add_option_group('MPFR Options')
    mpfr.add_option('--mpfr-dir',
                    help='Base directory where mpfr is installed')
    mpfr.add_option('--mpfr-incdir',
                    help='Directory where mpfr include files are installed')
    mpfr.add_option('--mpfr-libdir',
                    help='Directory where mpfr library files are installed')
    mpfr.add_option('--mpfr-libs',
                    help='Names of the mpfr libraries without prefix or suffix\n'
                         '(e.g. "mpfr")')
