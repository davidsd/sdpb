#! /usr/bin/env python
# encoding: utf-8

def configure(conf):
    def get_param(varname,default):
        return getattr(Options.options,varname,'')or default

    # Find GMPXX
    if conf.options.gmpxx_dir:
        if not conf.options.gmpxx_incdir:
            conf.options.gmpxx_incdir=conf.options.gmpxx_dir + "/include"
        if not conf.options.gmpxx_libdir:
            conf.options.gmpxx_libdir=conf.options.gmpxx_dir + "/lib"

    if conf.options.gmpxx_incdir:
        gmpxx_incdir=[conf.options.gmpxx_incdir]
    else:
        gmpxx_incdir=[]
    if conf.options.gmpxx_libdir:
        gmpxx_libdir=[conf.options.gmpxx_libdir]
    else:
        gmpxx_libdir=[]

    if conf.options.gmpxx_libs:
        gmpxx_libs=conf.options.gmpxx_libs.split()
    else:
        gmpxx_libs=['gmpxx','gmp']

    conf.check_cxx(msg="Checking for Gmpxx",
                   header_name='gmpxx.h',
                   includes=gmpxx_incdir,
                   uselib_store='gmpxx',
                   libpath=gmpxx_libdir,
                   rpath=gmpxx_libdir,
                   lib=gmpxx_libs)

def options(opt):
    gmpxx=opt.add_option_group('Gmpxx Options')
    gmpxx.add_option('--gmpxx-dir',
                   help='Base directory where gmpxx is installed')
    gmpxx.add_option('--gmpxx-incdir',
                   help='Directory where gmpxx include files are installed')
    gmpxx.add_option('--gmpxx-libdir',
                   help='Directory where gmpxx library files are installed')
    gmpxx.add_option('--gmpxx-libs',
                   help='Names of the gmpxx libraries without prefix or suffix\n'
                   '(e.g. "gmpxx")')
