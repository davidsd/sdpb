#! /usr/bin/env python
# encoding: utf-8

def configure(conf):
    # Find Boost
    import os
    if not conf.options.boost_incdir and not conf.options.boost_dir:
        for d in ['BOOST_INCLUDE','BOOST_INCLUDE_DIR','BOOST_INC_DIR']:
            env_dir=os.getenv(d)
            if env_dir:
                conf.to_log('Setting boost_incdir using environment variable: ' + d + '=' + env_dir)
                conf.options.boost_incdir=env_dir
                
    if not conf.options.boost_libdir and not conf.options.boost_dir:
        for d in ['BOOST_LIB','BOOST_LIB_DIR']:
            env_dir=os.getenv(d)
            if env_dir:
                conf.to_log('Setting boost_libdir using environment variable: ' + d + '=' + env_dir)
                conf.options.boost_libdir=env_dir

    if not conf.options.boost_dir:
        for d in ['BOOST_HOME','BOOST_DIR']:
            env_dir=os.getenv(d)
            if env_dir:
                conf.to_log('Setting boost_dir using environment variable: ' + d + '=' + env_dir)
                conf.options.boost_dir=env_dir
                
    if conf.options.boost_dir:
        if not conf.options.boost_incdir:
            conf.options.boost_incdir=conf.options.boost_dir + "/include"
        if not conf.options.boost_libdir:
            conf.options.boost_libdir=conf.options.boost_dir + "/lib"

    if conf.options.boost_incdir:
        boost_incdir=conf.options.boost_incdir.split()
    else:
        boost_incdir=[]
    if conf.options.boost_libdir:
        boost_libdir=conf.options.boost_libdir.split()
    else:
        boost_libdir=[]

    if conf.options.boost_libs:
        boost_libs=conf.options.boost_libs.split()
    else:
        boost_libs = ['boost_system', 'boost_date_time', 'boost_filesystem',
                      'boost_program_options', 'boost_iostreams', 'boost_serialization']

    boost_stacktrace_lib_found = any(name.startswith('boost_stacktrace') for name in boost_libs)
    # link to boost_stacktrace library instead of header-only compilation:
    boost_defines = ['BOOST_STACKTRACE_LINK'] if boost_stacktrace_lib_found else []

    conf.check_cxx(msg="Checking for Boost",
                   fragment="""#include <boost/iostreams/filter/gzip.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/device/file.hpp>
#include <boost/date_time/posix_time/posix_time.hpp>
#include <boost/program_options.hpp>
#include <boost/process.hpp>
#include <boost/serialization/serialization.hpp>
#include <boost/stacktrace.hpp>
int main()
{
boost::posix_time::second_clock::local_time();
boost::program_options::options_description();
boost::iostreams::file_sink("foo");
boost::iostreams::filtering_ostream();
boost::iostreams::gzip_compressor();
boost::process::ipstream pipe_stream;
boost::process::search_path("unzip");
boost::serialization::version_type version;
boost::stacktrace::stacktrace();
}
""",
                   includes=boost_incdir,
                   uselib_store='boost',
                   libpath=boost_libdir,
                   rpath=boost_libdir,
                   lib=boost_libs,
                   use=['cxx17'],
                   defines=boost_defines)

    # If boost_stacktrace library not defined by user, try linking to one of the libraries
    # listed in https://www.boost.org/doc/libs/1_84_0/doc/html/stacktrace/configuration_and_build.html
    # boost_stacktrace_backtrace and boost_stacktrace_addr2line can print source code location for each frame.
    # If all libraries fail to link, boost_stacktrace will be used as a header-only library
    # (which compiles longer and will not print source code location).
    if not boost_stacktrace_lib_found:
        for boost_stacktrace_lib in ['boost_stacktrace_backtrace', 'boost_stacktrace_addr2line',
                                     'boost_stacktrace_basic']:
            if conf.check_cxx(msg=f'Checking for {boost_stacktrace_lib}',
                              fragment="""
    #include <boost/stacktrace.hpp>
    int main()
    {
      boost::stacktrace::stacktrace();
    }
    """,
                              includes=boost_incdir,
                              uselib_store='boost',
                              libpath=boost_libdir,
                              rpath=boost_libdir,
                              lib=boost_libs + [boost_stacktrace_lib],
                              use=['cxx17'],
                              defines='BOOST_STACKTRACE_LINK',
                              mandatory=False):
                break



def options(opt):
    boost=opt.add_option_group('Boost Options')
    boost.add_option('--boost-dir',
                   help='Base directory where boost is installed')
    boost.add_option('--boost-incdir',
                   help='Directory where boost include files are installed')
    boost.add_option('--boost-libdir',
                   help='Directory where boost library files are installed')
    boost.add_option('--boost-libs',
                   help='Names of the boost libraries without prefix or suffix\n'
                        '(e.g. "boost_system boost_date_time boost_program_options")')
