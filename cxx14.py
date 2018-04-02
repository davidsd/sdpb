#! /usr/bin/env python
# encoding: utf-8

def configure(conf):
    def get_param(varname,default):
        return getattr(Options.options,varname,'')or default


    cxx14_fragment="#include <string>\nint main() {using namespace std::literals;\"hello\"s; }"
    flags=['','-std=c++14','-std=c++1y']
    if conf.options.cxx14_flag:
        flags=[conf.options.cxx14_flag]
    found_cxx14=False
    for flag in flags:
        try:
            conf.check_cxx(msg="Checking C++ 14 flag " + flag,
                           fragment=cxx14_fragment,
                           cxxflags=flag, linkflags=flag, uselib_store='cxx14')
        except conf.errors.ConfigurationError:
            continue
        else:
            found_cxx14=True
            break
    if not found_cxx14:
        conf.fatal('Could not find C++ 14 flag')

def options(opt):
    cxx14=opt.add_option_group('C++ 14 Options')
    cxx14.add_option('--cxx14-flag',
                     help='Flag to enable C++ 14')
