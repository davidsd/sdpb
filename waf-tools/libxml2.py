#! /usr/bin/env python
# encoding: utf-8

from check_config import check_config


def configure(conf):
    # Find Libxml2
    import os

    if not conf.options.libxml2_incdir:
        for d in ['LIBXML2_INCLUDE', 'LIBXML2_INCLUDE_DIR', 'LIBXML2_INC_DIR']:
            env_dir = os.getenv(d)
            if env_dir:
                conf.to_log('Setting libxml2_incdir using environment variable: ' + d + '=' + env_dir)
                conf.options.libxml2_incdir = env_dir

    if not conf.options.libxml2_libdir:
        for d in ['LIBXML2_LIB', 'LIBXML2_LIB_DIR']:
            env_dir = os.getenv(d)
            if env_dir:
                conf.to_log('Setting libxml2_libdir using environment variable: ' + d + '=' + env_dir)
                conf.options.libxml2_libdir = env_dir

    if conf.options.libxml2_dir:
        if not conf.options.libxml2_incdir:
            conf.options.libxml2_incdir = conf.options.libxml2_dir + "/include"
        if not conf.options.libxml2_libdir:
            conf.options.libxml2_libdir = conf.options.libxml2_dir + "/lib"

    if conf.options.libxml2_incdir:
        libxml2_incdir = conf.options.libxml2_incdir.split()
    else:
        libxml2_incdir = []
    if conf.options.libxml2_libdir:
        libxml2_libdir = conf.options.libxml2_libdir.split()
    else:
        libxml2_libdir = []

    if conf.options.libxml2_libs:
        libxml2_libs = conf.options.libxml2_libs.split()
    else:
        libxml2_libs = ['xml2']

    check_config(conf,
                 fragment='''
#include <libxml/parser.h>
int main()
{
  xmlInitParser();
}
      ''',
                 includes=libxml2_incdir,
                 uselib_store='libxml2',
                 libpath=libxml2_libdir,
                 lib=libxml2_libs,
                 packages=['libxml-2.0'])


def options(opt):
    libxml2 = opt.add_option_group('Libxml2 Options')
    libxml2.add_option('--libxml2-dir',
                       help='Base directory where libxml2 is installed')
    libxml2.add_option('--libxml2-incdir',
                       help='Directory where libxml2 include files are installed')
    libxml2.add_option('--libxml2-libdir',
                       help='Directory where libxml2 library files are installed')
    libxml2.add_option('--libxml2-libs',
                       help='Names of the libxml2 libraries without prefix or suffix\n'
                            '(e.g. "xml2")')
