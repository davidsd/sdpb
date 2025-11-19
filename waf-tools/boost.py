#! /usr/bin/env python
# encoding: utf-8

from check_config import check_config


def configure(conf):
    # Find Boost
    import os
    if not conf.options.boost_incdir and not conf.options.boost_dir:
        for d in ['BOOST_INCLUDE', 'BOOST_INCLUDE_DIR', 'BOOST_INC_DIR']:
            env_dir = os.getenv(d)
            if env_dir:
                conf.to_log('Setting boost_incdir using environment variable: ' + d + '=' + env_dir)
                conf.options.boost_incdir = env_dir

    if not conf.options.boost_libdir and not conf.options.boost_dir:
        for d in ['BOOST_LIB', 'BOOST_LIB_DIR']:
            env_dir = os.getenv(d)
            if env_dir:
                conf.to_log('Setting boost_libdir using environment variable: ' + d + '=' + env_dir)
                conf.options.boost_libdir = env_dir

    if not conf.options.boost_dir:
        for d in ['BOOST_HOME', 'BOOST_DIR', 'BOOST_ROOT']:
            env_dir = os.getenv(d)
            if env_dir:
                conf.to_log('Setting boost_dir using environment variable: ' + d + '=' + env_dir)
                conf.options.boost_dir = env_dir

    if conf.options.boost_dir:
        if not conf.options.boost_incdir:
            conf.options.boost_incdir = conf.options.boost_dir + "/include"
        if not conf.options.boost_libdir:
            conf.options.boost_libdir = conf.options.boost_dir + "/lib"

    if conf.options.boost_incdir:
        boost_incdir = conf.options.boost_incdir.split()
    else:
        boost_incdir = []
    if conf.options.boost_libdir:
        boost_libdir = conf.options.boost_libdir.split()
    else:
        boost_libdir = []

    if conf.options.boost_libs:
        boost_libs = conf.options.boost_libs.split()
    else:
        boost_libs = ['boost_date_time', 'boost_filesystem',
                      'boost_program_options', 'boost_iostreams', 'boost_serialization']

    # Boost.System library is header-only since 1.69
    # libboost_system stub has been removed since 1.89
    boost_system_lib = 'boost_system'
    boost_system_lib_found = any(name == boost_system_lib for name in boost_libs)

    # Boost.Process v2
    boost_process_lib = 'boost_process'
    boost_process_lib_found = any(name == boost_process_lib for name in boost_libs)

    boost_stacktrace_lib_found = any(name.startswith('boost_stacktrace') for name in boost_libs)
    # link to boost_stacktrace library instead of header-only compilation:
    boost_defines = ['BOOST_STACKTRACE_LINK'] if boost_stacktrace_lib_found else []

    conf.start_msg('Checking for Boost')

    # Link to boost_system, if necessary
    if not boost_system_lib_found:
        for boost_system_libs in [[], [boost_system_lib]]:
            if check_config(conf,
                            msg=f'  Checking for Boost.System, boost_system_libs={boost_system_libs}',
                            fragment="""
#include <boost/system/error_code.hpp>
namespace sys = boost::system;
int main()
{
  sys::error_code ec;
  const auto& cat = sys::generic_category();
}
""",
                            includes=boost_incdir,
                            uselib_store='boost',
                            libpath=boost_libdir,
                            lib=boost_system_libs,
                            use=['cxx17'],
                            mandatory=False):
                boost_system_lib_found = True
                boost_libs += boost_system_libs
                break

    if not boost_system_lib_found:
        conf.fatal('Could not find Boost.System')

    # Check other Boost libraries
    check_config(conf,
                 fragment="""#include <boost/iostreams/filter/gzip.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/device/file.hpp>
#include <boost/date_time/posix_time/posix_time.hpp>
#include <boost/program_options.hpp>
#include <boost/serialization/serialization.hpp>
#include <boost/stacktrace.hpp>
int main()
{
boost::posix_time::second_clock::local_time();
boost::program_options::options_description();
boost::iostreams::file_sink("foo");
boost::iostreams::filtering_ostream();
boost::iostreams::gzip_compressor();
boost::serialization::version_type version;
boost::stacktrace::stacktrace();
}
""",
                 includes=boost_incdir,
                 uselib_store='boost',
                 libpath=boost_libdir,
                 lib=boost_libs,
                 use=['cxx17'],
                 defines=boost_defines)

    # Link to boost_process, if necessary
    if not boost_process_lib_found:
        for boost_process_libs in [[], [boost_process_lib]]:
            if check_config(conf,
                            msg=f'  Checking for Boost.Process, boost_process_libs={boost_process_libs}',
                            fragment="""
#include <boost/version.hpp>

// There are two versions for Boost.Process, V1 and V2.
// V2 was introduced in Boost 1.80 and became default in Boost 1.89.
// Prior to Boost 1.86, V2 failed to compile on Alpine Linux (musl-libc),
// see https://github.com/boostorg/process/pull/376
// NB: The code for choosing Boost.Process version should match test/src/integration_tests/util/process.cxx
#if BOOST_VERSION < 108600
#define USE_BOOST_PROCESS_V1
#include <boost/process.hpp>
namespace bp = boost::process;

#elif BOOST_VERSION < 108900
#include <boost/process/v2.hpp>
namespace bp = boost::process::v2;

#else
#include <boost/process.hpp>
namespace bp = boost::process;

#endif

#include <filesystem>

std::filesystem::path find_executable(const std::string &filename)
{
#ifdef USE_BOOST_PROCESS_V1
  return bp::search_path(filename).string();
#else
  return bp::environment::find_executable(filename).string();
#endif
}

int main()
{
  auto ls = find_executable("ls");
}
""",
                            includes=boost_incdir,
                            uselib_store='boost',
                            libpath=boost_libdir,
                            lib=boost_libs + boost_process_libs,
                            use=['cxx17'],
                            mandatory=False):
                boost_process_lib_found = True
                boost_libs += boost_process_libs
                break

    if not boost_process_lib_found:
        conf.fatal('Could not find Boost.Process')

    # If boost_stacktrace library not defined by user, try linking to one of the libraries
    # listed in https://www.boost.org/doc/libs/1_84_0/doc/html/stacktrace/configuration_and_build.html
    # boost_stacktrace_backtrace and boost_stacktrace_addr2line can print source code location for each frame.
    # If all libraries fail to link, boost_stacktrace will be used as a header-only library
    # (which compiles longer and will not print source code location).
    if not boost_stacktrace_lib_found:
        for boost_stacktrace_lib in ['boost_stacktrace_backtrace', 'boost_stacktrace_addr2line',
                                     'boost_stacktrace_basic']:
            if check_config(conf,
                            msg='  Checking for ' + boost_stacktrace_lib,
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
                            lib=boost_libs + [boost_stacktrace_lib],
                            use=['cxx17'],
                            defines='BOOST_STACKTRACE_LINK',
                            mandatory=False):
                break

    conf.end_msg(True)


def options(opt):
    boost = opt.add_option_group('Boost Options')
    boost.add_option('--boost-dir',
                     help='Base directory where boost is installed')
    boost.add_option('--boost-incdir',
                     help='Directory where boost include files are installed')
    boost.add_option('--boost-libdir',
                     help='Directory where boost library files are installed')
    boost.add_option('--boost-libs',
                     help='Names of the boost libraries without prefix or suffix\n'
                          '(e.g. "boost_system boost_date_time boost_program_options")')
