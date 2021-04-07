#! /usr/bin/env python
# encoding: utf-8

def configure(conf):
    # Find Libarchive_Cpp
    import os
    if not conf.options.libarchive_cpp_incdir:
        for d in ['LIBARCHIVE_CPP_INCLUDE','LIBARCHIVE_CPP_INCLUDE_DIR','LIBARCHIVE_CPP_INC_DIR']:
            env_dir=os.getenv(d)
            if env_dir:
                conf.to_log('Setting libarchive_cpp_incdir using environment variable: ' + d + '=' + env_dir)
                conf.options.libarchive_cpp_incdir=env_dir
                
    if not conf.options.libarchive_cpp_libdir:
        for d in ['LIBARCHIVE_CPP_LIB','LIBARCHIVE_CPP_LIB_DIR']:
            env_dir=os.getenv(d)
            if env_dir:
                conf.to_log('Setting libarchive_cpp_libdir using environment variable: ' + d + '=' + env_dir)
                conf.options.libarchive_cpp_libdir=env_dir

    if conf.options.libarchive_cpp_dir:
        if not conf.options.libarchive_cpp_incdir:
            conf.options.libarchive_cpp_incdir=conf.options.libarchive_cpp_dir + "/include"
        if not conf.options.libarchive_cpp_libdir:
            conf.options.libarchive_cpp_libdir=conf.options.libarchive_cpp_dir + "/lib"

    if conf.options.libarchive_cpp_incdir:
        libarchive_cpp_incdir=conf.options.libarchive_cpp_incdir.split()
    else:
        libarchive_cpp_incdir=[]
    if conf.options.libarchive_cpp_libdir:
        libarchive_cpp_libdir=conf.options.libarchive_cpp_libdir.split()
    else:
        libarchive_cpp_libdir=[]

    if conf.options.libarchive_cpp_libs:
        libarchive_cpp_libs=conf.options.libarchive_cpp_libs.split()
    else:
        libarchive_cpp_libs=['archive','archive_cpp_wrapper']

    if not conf.check_cxx(msg="Checking for libarchive_cpp",
                          fragment='''
#include <archive_reader.hpp>
#include <archive_exception.hpp>
#include <fstream>
int main()
{
  std::fstream fs;
  ns_archive::reader(
    ns_archive::reader::make_reader<ns_archive::ns_reader::format::_ALL, ns_archive::ns_reader::filter::_ALL>(
      fs, 10240));
  ;
}
''',
                          includes=libarchive_cpp_incdir,
                          uselib_store='libarchive_cpp',
                          libpath=libarchive_cpp_libdir,
                          rpath=libarchive_cpp_libdir,
                          lib=libarchive_cpp_libs,
                          use=['cxx14']):
        conf.fatal("Could not find libarchive_cpp")


def options(opt):
    libarchive_cpp=opt.add_option_group('Libarchive_Cpp Options')
    libarchive_cpp.add_option('--libarchive_cpp-dir',
                   help='Base directory where libarchive_cpp is installed')
    libarchive_cpp.add_option('--libarchive_cpp-incdir',
                   help='Directory where libarchive_cpp include files are installed')
    libarchive_cpp.add_option('--libarchive_cpp-libdir',
                   help='Directory where libarchive_cpp library files are installed')
    libarchive_cpp.add_option('--libarchive_cpp-libs',
                   help='Names of the libarchive and libarchive_cpp libraries without prefix or suffix\n'
                   '(e.g. "archive archive_cpp")')
