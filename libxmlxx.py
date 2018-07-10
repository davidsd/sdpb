#! /usr/bin/env python
# encoding: utf-8

def configure(conf):
    conf.load('compiler_cxx cxx14')
    
    # Find Libxmlxx
    if conf.options.libxmlxx_dir:
        if not conf.options.libxmlxx_incdir:
            conf.options.libxmlxx_incdir=conf.options.libxmlxx_dir + "/include"
        if not conf.options.libxmlxx_libdir:
            conf.options.libxmlxx_libdir=conf.options.libxmlxx_dir + "/lib"

    if conf.options.libxmlxx_incdir:
        libxmlxx_incdir=conf.options.libxmlxx_incdir.split()
    else:
        libxmlxx_incdir=[]
    if conf.options.libxmlxx_libdir:
        libxmlxx_libdir=conf.options.libxmlxx_libdir.split()
    else:
        libxmlxx_libdir=[]

    if conf.options.libxmlxx_libs:
        libxmlxx_libs=conf.options.libxmlxx_libs.split()
    else:
        libxmlxx_libs=['xml++-2.6', 'xml2', 'glibmm-2.4', 'gobject-2.0',
                       'glib-2.0', 'sigc-2.0']

    if not conf.check_cxx(msg="Checking for libxml++",
                          header_name='libxml++/libxml++.h',
                          includes=libxmlxx_incdir,
                          uselib_store='libxmlxx',
                          libpath=libxmlxx_libdir,
                          rpath=libxmlxx_libdir,
                          lib=libxmlxx_libs,
                          mandatory=False) and not conf.check_cfg(path='pkg-config', args='--libs --cflags',
                                                                  package='libxml++-2.6', uselib_store='libxmlxx',
                                                                  mandatory=False) and not conf.check_cfg(path='pkg-config', args='--libs --cflags',
		                                                                                          package='libxml++-', uselib_store='libxmlxx',
                                                                                                          mandatory=True):
        conf.fatal("Could not find libxml++.")


def options(opt):
    libxmlxx=opt.add_option_group('Libxmlxx Options')
    libxmlxx.add_option('--libxmlxx-dir',
                   help='Base directory where libxmlxx is installed')
    libxmlxx.add_option('--libxmlxx-incdir',
                   help='Directory where libxmlxx include files are installed')
    libxmlxx.add_option('--libxmlxx-libdir',
                   help='Directory where libxmlxx library files are installed')
    libxmlxx.add_option('--libxmlxx-libs',
                   help='Names of the libxmlxx libraries without prefix or suffix\n'
                   '(e.g. "xml++-2.6 xml2 glibmm-2.4 gobject-2.0 glib-2.0 sigc-2.0")')
