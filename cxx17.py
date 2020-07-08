#! /usr/bin/env python
# encoding: utf-8

def configure(conf):
    def get_param(varname,default):
        return getattr(Options.options,varname,'')or default


    cxx17_fragment="#include <string>\n#include <optional>\nint main() {std::optional<int> i;\nusing namespace std::literals;\"hello\"s; }"
    flags=['','-std=c++17','-std=c++1z']
    if conf.options.cxx17_flag:
        flags=[conf.options.cxx17_flag]
    found_cxx17=False
    for flag in flags:
        try:
            conf.check_cxx(msg="Checking C++ 17 flag " + flag,
                           fragment=cxx17_fragment,
                           cxxflags=flag, linkflags=flag, uselib_store='cxx17')
        except conf.errors.ConfigurationError:
            continue
        else:
            found_cxx17=True
            break
    if not found_cxx17:
        conf.fatal('Could not find C++ 17 flag')

def options(opt):
    cxx17=opt.add_option_group('C++ 17 Options')
    cxx17.add_option('--cxx17-flag',
                     help='Flag to enable C++ 17')
