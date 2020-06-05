#! /usr/bin/env python
# encoding: utf-8

def configure(conf):
    def get_param(varname,default):
        return getattr(Options.options,varname,'')or default

    import os
    # Find Rapidjson
    if conf.options.rapidjson_dir:
        if not conf.options.rapidjson_incdir:
            conf.options.rapidjson_incdir=conf.options.rapidjson_dir + "/include"

    if conf.options.rapidjson_incdir:
        rapidjson_incdir=conf.options.rapidjson_incdir.split()
    else:
        rapidjson_incdir=[]

    conf.check_cxx(msg="Checking for rapidjson",
                   fragment="#include <rapidjson/reader.h>\nint main() {rapidjson::Reader();}\n",
                   includes=rapidjson_incdir,
                   uselib_store='rapidjson',
                   use=['cxx14'])

def options(opt):
    rapidjson=opt.add_option_group('Rapidjson Options')
    rapidjson.add_option('--rapidjson-dir',
                   help='Base directory where rapidjson is installed')
    rapidjson.add_option('--rapidjson-incdir',
                   help='Directory where rapidjson include files are installed')
