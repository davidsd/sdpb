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
        boost_libs=['boost_system', 'boost_filesystem', 'boost_date_time',
                    'boost_program_options', 'boost_iostreams']

    conf.check_cxx(msg="Checking for Boost",
                   fragment="""#include <boost/iostreams/filter/gzip.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/device/file.hpp>
#include <boost/filesystem.hpp>
#include <boost/date_time/posix_time/posix_time.hpp>
#include <boost/program_options.hpp>
int main()
{
boost::posix_time::second_clock::local_time();
boost::filesystem::path();
boost::program_options::options_description();
boost::iostreams::file_sink("foo");
boost::iostreams::filtering_ostream();
boost::iostreams::gzip_compressor();
}
""",
                   includes=boost_incdir,
                   uselib_store='boost',
                   libpath=boost_libdir,
                   rpath=boost_libdir,
                   lib=boost_libs,
                   use=['cxx14'])



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
                   '(e.g. "boost_system boost_filesystem boost_date_time boost_program_options")')
