#! /usr/bin/env python
# encoding: utf-8

import os

def configure(conf):
    def get_param(varname,default):
        return getattr(Options.options,varname,'')or default

    if not conf.options.apfp_incdir:
        for d in ['APFP_INCLUDE','APFP_INCLUDE_DIR','APFP_INC_DIR']:
            env_dir=os.getenv(d)
            if env_dir:
                conf.to_log('Setting apfp_incdir using environment variable: ' + d + '=' + env_dir)
                conf.options.apfp_incdir=env_dir
                
    if not conf.options.apfp_libdir:
        for d in ['APFP_LIB','APFP_LIB_DIR']:
            env_dir=os.getenv(d)
            if env_dir:
                conf.to_log('Setting apfp_libdir using environment variable: ' + d + '=' + env_dir)
                conf.options.apfp_libdir=env_dir

    if not conf.options.apfp_dir:
        env_dir=os.getenv('APFP_DIR')
        if env_dir:
            conf.to_log('Setting apfp_dir using environment variable: ' + d + '=' + env_dir)
            conf.options.apfp_dir=env_dir
                
    # Find APFP
    if conf.options.apfp_dir:
        conf.to_log('Using apfp_dir: ' + conf.options.apfp_dir)
        if not conf.options.apfp_incdir:
            conf.to_log('Setting apfp_incdir using apfp_dir: ' + conf.options.apfp_dir)
            conf.options.apfp_incdir=conf.options.apfp_dir + "/include"
        if not conf.options.apfp_libdir:
            conf.to_log('Setting apfp_libdir using apfp_dir: ' + conf.options.apfp_dir)
            conf.options.apfp_libdir=conf.options.apfp_dir + "/lib"

    if conf.options.apfp_incdir:
        conf.to_log('Using apfp_incdir: ' + conf.options.apfp_incdir)
        apfp_incdir=conf.options.apfp_incdir.split()
    else:
        apfp_incdir=[]
    if conf.options.apfp_libdir:
        conf.to_log('Using apfp_libdir: ' + conf.options.apfp_libdir)
        apfp_libdir=conf.options.apfp_libdir.split()
    else:
        apfp_libdir=[]

    if conf.options.apfp_libs:
        conf.to_log('Using apfp_libs: ' + conf.options.apfp_libs)
        apfp_libs=conf.options.apfp_libs.split()
    else:
        apfp_libs=['apfp', 'dl', 'stdc++fs', 'OpenCL', 'xilinxopencl', 'xrt_core', 'pthread']

    conf.check_cxx(msg="Checking for APFP",
                   header_name='apfp/config.h',
                   includes=apfp_incdir,
                   uselib_store='apfp',
                   libpath=apfp_libdir,
                   rpath=apfp_libdir,
                   lib=apfp_libs,
                   use=['cxx14', 'gmpxx', 'mpfr'])

def options(opt):
    apfp=opt.add_option_group('APFP Options')
    apfp.add_option('--apfp-dir',
                   help='Base directory where apfp is installed')
    apfp.add_option('--apfp-incdir',
                   help='Directory where apfp include files are installed')
    apfp.add_option('--apfp-libdir',
                   help='Directory where apfp library files are installed')
    apfp.add_option('--apfp-libs',
                   help='Names of the apfp libraries without prefix or suffix\n'
                   '(e.g. "apfp")')
