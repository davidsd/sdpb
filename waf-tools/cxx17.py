#! /usr/bin/env python
# encoding: utf-8

def configure(conf):
    def get_param(varname,default):
        return getattr(Options.options,varname,'')or default

    cxx17_fragment='''
    #include <string>
    #include <optional>
    #include <filesystem>
    int main() {
        std::optional<int> i;
        using namespace std::literals;
        "hello\"s;
        std::filesystem::path path(".");
    }
    '''
    flag='-std=c++17'
    try:
        conf.check_cxx(msg="Checking for C++17 features",
                        fragment=cxx17_fragment,
                        cxxflags=flag, linkflags=flag, uselib_store='cxx17')
    except conf.errors.ConfigurationError as e:
        conf.fatal('Failed to check C++ 17 features')
