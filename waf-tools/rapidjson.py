#! /usr/bin/env python
# encoding: utf-8

from check_config import check_config


def configure(conf):
    # Find Rapidjson
    if conf.options.rapidjson_dir:
        if not conf.options.rapidjson_incdir:
            conf.options.rapidjson_incdir = conf.options.rapidjson_dir + "/include"

    if conf.options.rapidjson_incdir:
        rapidjson_incdir = conf.options.rapidjson_incdir.split()
    else:
        rapidjson_incdir = []

    check_config(conf,
                 fragment="#include <rapidjson/reader.h>\nint main() {rapidjson::Reader();}\n",
                 includes=rapidjson_incdir,
                 uselib_store='rapidjson',
                 libpath=[],
                 lib=[],
                 use=['cxx17'])


def options(opt):
    rapidjson = opt.add_option_group('Rapidjson Options')
    rapidjson.add_option('--rapidjson-dir',
                         help='Base directory where rapidjson is installed')
    rapidjson.add_option('--rapidjson-incdir',
                         help='Directory where rapidjson include files are installed')
