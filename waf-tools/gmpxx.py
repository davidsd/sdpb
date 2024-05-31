#! /usr/bin/env python
# encoding: utf-8

from check_config import check_config


def configure(conf):
    import os

    if not conf.options.gmpxx_incdir:
        for d in ['GMP_INCLUDE', 'GMPXX_INCLUDE', 'GMP_INCLUDE_DIR',
                  'GMPXX_INCLUDE_DIR', 'GMP_INC_DIR', 'GMPXX_INC_DIR']:
            env_dir = os.getenv(d)
            if env_dir:
                conf.to_log('Setting gmpxx_incdir using environment variable: ' + d + '=' + env_dir)
                conf.options.gmpxx_incdir = env_dir

    if not conf.options.gmpxx_libdir:
        for d in ['GMP_LIB', 'GMPXX_LIB', 'GMP_LIB_DIR', 'GMPXX_LIB_DIR']:
            env_dir = os.getenv(d)
            if env_dir:
                conf.to_log('Setting gmpxx_libdir using environment variable: ' + d + '=' + env_dir)
                conf.options.gmpxx_libdir = env_dir

    if not conf.options.gmpxx_dir:
        for d in ['GMP_DIR', 'GMPXX_DIR']:
            env_dir = os.getenv(d)
            if env_dir:
                conf.to_log('Setting gmpxx_dir using environment variable: ' + d + '=' + env_dir)
                conf.options.gmpxx_dir = env_dir

    # Find GMPXX
    if conf.options.gmpxx_dir:
        conf.to_log('Using gmpxx_dir: ' + conf.options.gmpxx_dir)
        if not conf.options.gmpxx_incdir:
            conf.to_log('Setting gmpxx_incdir using gmpxx_dir: ' + conf.options.gmpxx_dir)
            conf.options.gmpxx_incdir = conf.options.gmpxx_dir + "/include"
        if not conf.options.gmpxx_libdir:
            conf.to_log('Setting gmpxx_libdir using gmpxx_dir: ' + conf.options.gmpxx_dir)
            conf.options.gmpxx_libdir = conf.options.gmpxx_dir + "/lib"

    if conf.options.gmpxx_incdir:
        conf.to_log('Using gmpxx_incdir: ' + conf.options.gmpxx_incdir)
        gmpxx_incdir = conf.options.gmpxx_incdir.split()
    else:
        gmpxx_incdir = []
    if conf.options.gmpxx_libdir:
        conf.to_log('Using gmpxx_libdir: ' + conf.options.gmpxx_libdir)
        gmpxx_libdir = conf.options.gmpxx_libdir.split()
    else:
        gmpxx_libdir = []

    if conf.options.gmpxx_libs:
        conf.to_log('Using gmpxx_libs: ' + conf.options.gmpxx_libs)
        gmpxx_libs = conf.options.gmpxx_libs.split()
    else:
        gmpxx_libs = ['gmpxx', 'gmp']

    check_config(conf,
                 fragment='''
#include <gmpxx.h>
int main()
{
    mpf_t qf;
    mpf_init(qf);
}
                 ''',
                 includes=gmpxx_incdir,
                 uselib_store='gmpxx',
                 libpath=gmpxx_libdir,
                 lib=gmpxx_libs,
                 use=['cxx17'])


def options(opt):
    gmpxx = opt.add_option_group('GMPXX Options')
    gmpxx.add_option('--gmpxx-dir',
                     help='Base directory where gmpxx is installed')
    gmpxx.add_option('--gmpxx-incdir',
                     help='Directory where gmpxx include files are installed')
    gmpxx.add_option('--gmpxx-libdir',
                     help='Directory where gmpxx library files are installed')
    gmpxx.add_option('--gmpxx-libs',
                     help='Names of the gmpxx libraries without prefix or suffix\n'
                          '(e.g. "gmpxx")')
