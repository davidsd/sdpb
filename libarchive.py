#! /usr/bin/env python
# encoding: utf-8

def configure(conf):
    # Find Libarchive
    import os
    if not conf.options.libarchive_incdir:
        for d in ['LIBARCHIVE_INCLUDE','LIBARCHIVE_INCLUDE_DIR','LIBARCHIVE_INC_DIR']:
            env_dir=os.getenv(d)
            if env_dir:
                conf.to_log('Setting libarchive_incdir using environment variable: ' + d + '=' + env_dir)
                conf.options.libarchive_incdir=env_dir
                
    if not conf.options.libarchive_libdir:
        for d in ['LIBARCHIVE_LIB','LIBARCHIVE_LIB_DIR']:
            env_dir=os.getenv(d)
            if env_dir:
                conf.to_log('Setting libarchive_libdir using environment variable: ' + d + '=' + env_dir)
                conf.options.libarchive_libdir=env_dir

    if conf.options.libarchive_dir:
        if not conf.options.libarchive_incdir:
            conf.options.libarchive_incdir=conf.options.libarchive_dir + "/include"
        if not conf.options.libarchive_libdir:
            conf.options.libarchive_libdir=conf.options.libarchive_dir + "/lib"

    if conf.options.libarchive_incdir:
        libarchive_incdir=conf.options.libarchive_incdir.split()
    else:
        libarchive_incdir=[]
    if conf.options.libarchive_libdir:
        libarchive_libdir=conf.options.libarchive_libdir.split()
    else:
        libarchive_libdir=[]

    if conf.options.libarchive_libs:
        libarchive_libs=conf.options.libarchive_libs.split()
    else:
        libarchive_libs=['archive']

    if not conf.check_cxx(msg="Checking for libarchive",
                          fragment='''
#include <archive.h>
int main()
{
  archive *a(archive_read_new());
  archive_read_free(a);
}
''',
                          includes=libarchive_incdir,
                          uselib_store='libarchive',
                          libpath=libarchive_libdir,
                          rpath=libarchive_libdir,
                          lib=libarchive_libs,
                          use=['cxx14']):
        conf.fatal("Could not find libarchive")


def options(opt):
    libarchive=opt.add_option_group('Libarchive Options')
    libarchive.add_option('--libarchive-dir',
                   help='Base directory where libarchive is installed')
    libarchive.add_option('--libarchive-incdir',
                   help='Directory where libarchive include files are installed')
    libarchive.add_option('--libarchive-libdir',
                   help='Directory where libarchive library files are installed')
    libarchive.add_option('--libarchive-libs',
                   help='Names of the libarchive and libarchive libraries without prefix or suffix\n'
                   '(e.g. "archive")')
