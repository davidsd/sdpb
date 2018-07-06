#! /usr/bin/env python
# encoding: utf-8

import os

def configure(conf):
    def get_param(varname,default):
        return getattr(Options.options,varname,'')or default

    if not conf.options.gmpxx_incdir:
        for d in ['GMP_INCLUDE','GMPXX_INCLUDE','GMP_INCLUDE_DIR',
                  'GMPXX_INCLUDE_DIR','GMP_INC_DIR','GMPXX_INC_DIR']:
            env_dir=os.getenv(d)
            if env_dir:
                conf.options.gmpxx_incdir=env_dir
                
    if not conf.options.gmpxx_libdir:
        for d in ['GMP_LIB','GMPXX_LIB','GMP_LIB_DIR','GMPXX_LIB_DIR']:
            env_dir=os.getenv(d)
            if env_dir:
                conf.options.gmpxx_incdir=env_dir

    if not conf.options.gmpxx_dir:
        for d in ['GMP_DIR','GMPXX_DIR']:
            env_dir=os.getenv(d)
            if env_dir:
                conf.options.gmpxx_incdir=env_dir
                
    # Find GMPXX
    if conf.options.gmpxx_dir:
        if not conf.options.gmpxx_incdir:
            conf.options.gmpxx_incdir=conf.options.gmpxx_dir + "/include"
        if not conf.options.gmpxx_libdir:
            conf.options.gmpxx_libdir=conf.options.gmpxx_dir + "/lib"

    if conf.options.gmpxx_incdir:
        gmpxx_incdir=conf.options.gmpxx_incdir.split()
    else:
        gmpxx_incdir=[]
    if conf.options.gmpxx_libdir:
        gmpxx_libdir=conf.options.gmpxx_libdir.split()
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
