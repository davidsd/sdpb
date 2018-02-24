#! /usr/bin/env python
# encoding: utf-8

def configure(conf):
    def get_param(varname,default):
        return getattr(Options.options,varname,'')or default

    # Find TINYXML2
    if conf.options.tinyxml2_dir:
        if not conf.options.tinyxml2_incdir:
            conf.options.tinyxml2_incdir=conf.options.tinyxml2_dir + "/include"
        if not conf.options.tinyxml2_libdir:
            conf.options.tinyxml2_libdir=conf.options.tinyxml2_dir + "/lib"

    if conf.options.tinyxml2_incdir:
        tinyxml2_incdir=[conf.options.tinyxml2_incdir]
    else:
        tinyxml2_incdir=[]
    if conf.options.tinyxml2_libdir:
        tinyxml2_libdir=[conf.options.tinyxml2_libdir]
    else:
        tinyxml2_libdir=[]

    if conf.options.tinyxml2_libs:
        tinyxml2_libs=conf.options.tinyxml2_libs.split()
    else:
        tinyxml2_libs=['tinyxml2']

    conf.check_cxx(msg="Checking for Tinyxml2",
                   header_name='tinyxml2.h',
                   includes=tinyxml2_incdir,
                   uselib_store='tinyxml2',
                   libpath=tinyxml2_libdir,
                   rpath=tinyxml2_libdir,
                   lib=tinyxml2_libs)

def options(opt):
    tinyxml2=opt.add_option_group('Tinyxml2 Options')
    tinyxml2.add_option('--tinyxml2-dir',
                   help='Base directory where tinyxml2 is installed')
    tinyxml2.add_option('--tinyxml2-incdir',
                   help='Directory where tinyxml2 include files are installed')
    tinyxml2.add_option('--tinyxml2-libdir',
                   help='Directory where tinyxml2 library files are installed')
    tinyxml2.add_option('--tinyxml2-libs',
                   help='Names of the tinyxml2 libraries without prefix or suffix\n'
                   '(e.g. "tinyxml2")')
