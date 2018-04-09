#! /usr/bin/env python

def configure(conf):
    def get_param(varname,default):
        return getattr(Options.options,varname,'')or default

    # Find CUDA
    if(conf.options.cuda):
        if conf.options.cuda_dir:
            if not conf.options.cuda_incdir:
                conf.options.cuda_incdir=conf.options.cuda_dir + "/include"
            if not conf.options.cuda_libdir:
                conf.options.cuda_libdir=conf.options.cuda_dir + "/lib"

        if conf.options.cuda_incdir:
            cuda_incdir=[conf.options.cuda_incdir]
        else:
            cuda_incdir=[]
        if conf.options.cuda_libdir:
            cuda_libdir=[conf.options.cuda_libdir]
        else:
            cuda_libdir=[]

        if conf.options.cuda_libs:
            cuda_libs=conf.options.cuda_libs.split()
        else:
            cuda_libs=['cudart','curand','cublas']

        cuda_incdir=[]
        cuda_libdir=[]
        cuda_cxxflags=['-D__SDPB_CUDA__']
        cuda_linkflags=[]
        cuda_libs=['cudart','cublas','curand']

        conf.check_cxx(msg="Checking for CUDA",
                       header_name='cuda.h',
                       includes=cuda_incdir,
                       cxxflags=cuda_cxxflags,
                       linkflags=cuda_linkflags,
                       uselib_store='cuda',
                       libpath=cuda_libdir,
                       rpath=cuda_libdir,
                       lib=cuda_libs)

def options(opt):
    cuda=opt.add_option_group('CUDA Options')
    cuda.add_option('--cuda',default=False,action='store_true',
                   help='Use CUDA/CUBLAS where possible')
    cuda.add_option('--cuda-dir',
                   help='Base directory where CUDA is installed')
    cuda.add_option('--cuda-incdir',
                   help='Directory where CUDA include files are installed')
    cuda.add_option('--cuda-cxxflags',
                   help='Additional flags when compiling (e.g. -Dcuda_ILP64 -m64)')
    cuda.add_option('--cuda-libdir',
                   help='Directory where CUDA library files are installed')
    cuda.add_option('--cuda-libs',
                   help='Names of the CUDA libraries without prefix or suffix\n'
                   '(e.g. "cudart")')
    cuda.add_option('--cuda-linkflags',
                   help='Additional flags when linking (e.g. -Wl,--no-as-needed)')
