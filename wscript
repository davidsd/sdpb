import os, subprocess
from waflib import TaskGen

# Include and library paths are searched in this order
external_packages = ['cxx17', 'gmpxx', 'mpfr', 'boost', 'elemental', 'libxml2', 'rapidjson', 'libarchive', 'flint',
                     'cblas']


@TaskGen.feature('cxxprogram')
@TaskGen.after_method('apply_link')
def relink_when_gpu_library_changes(self):
    # the archive is passed in LDFLAGS, which waf does not track as a dependency
    if self.env.SDPB_GPU and getattr(self, 'gpu_library', False):
        self.link_task.dep_nodes.append(self.bld.bldnode.find_or_declare('libsdpb_gpu.a'))


def options(opt):
    opt.load(['compiler_cxx', 'gnu_dirs'])
    opt.add_option('--enable-gpu', action='store_true', default=False,
                   help='Build sdpb with the GPU backend (src/gpu, CUDA): the --gpu option of sdpb')
    opt.add_option('--nvcc', action='store', default='nvcc',
                   help='CUDA compiler, with options, e.g. "nvcc --compiler-bindir /path/to/gcc13/bin"')
    opt.add_option('--gpu-arch', action='store', default='90', help='CUDA architectures, comma separated (90: H100/H200)')
    opt.add_option('--gpu-nccl-dir', action='store', default='', help='NCCL headers for multi-GPU support (library loaded at run time)')
    opt.add_option('--cuda-dir', action='store', default=os.environ.get('CUDA_HOME', ''),
                   help='CUDA toolkit prefix (lib64/ has libcudart_static.a and libcublasLt)')
    opt.add_option('--gpu-gmp-include', action='store', default='',
                   help='Directory of gmp.h for the CUDA compiler')
    opt.load(
        ['cxx17', 'boost', 'gmpxx', 'mpfr', 'mpsolve', 'elemental', 'libxml2', 'rapidjson', 'libarchive', 'cblas',
         'flint'],
        tooldir='./waf-tools')


def configure(conf):
    if not 'CXX' in os.environ or os.environ['CXX'] == 'g++' or os.environ['CXX'] == 'icpc':
        conf.environ['CXX'] = 'mpicxx'

    conf.load(['compiler_cxx', 'gnu_dirs', 'cxx17', 'boost', 'gmpxx', 'mpfr', 'mpsolve',
               'elemental', 'libxml2', 'rapidjson', 'libarchive', 'cblas', 'flint'])
    conf.load('clang_compilation_database', tooldir='./waf-tools')
    conf.load('check_package_versions', tooldir='./waf-tools')
    # Some targets also use mpsolve, after all other packages
    conf.check_package_versions(external_packages + ['mpsolve'])

    conf.env.git_version = subprocess.check_output('git describe --tags --always --dirty', universal_newlines=True,
                                                   shell=True).rstrip()

    # GPU backend, off by default; without it the CPU build is unchanged
    conf.env.SDPB_GPU = bool(conf.options.enable_gpu)
    if conf.env.SDPB_GPU:
        conf.env.NVCC = conf.options.nvcc
        conf.env.GPU_ARCH = conf.options.gpu_arch.split(',')
        conf.env.CUDA_DIR = conf.options.cuda_dir
        conf.env.GPU_GMP_INCLUDE = conf.options.gpu_gmp_include
        conf.env.GPU_NCCL_DIR = conf.options.gpu_nccl_dir


def write_version_header(task):
    task.outputs[0].write('#pragma once\n#define SDPB_VERSION_STRING "%s"\n' % task.env.git_version)


def build(bld):
    # Generate build/sdpb_version.hxx before compiling anything
    bld(rule=write_version_header, target='sdpb_version.hxx', vars=['git_version'])
    bld.add_group()

    default_flags = ['-Wall', '-Wextra', '-Werror=return-type', '-O3', '-g']
    default_defines = ['OMPI_SKIP_MPICXX']
    gpu_link = {}
    if bld.env.SDPB_GPU:
        # The device library (src/gpu) is compiled by nvcc and linked into sdpb.
        default_defines = default_defines + ['SDPB_GPU']
        gpu_sources = ['src/gpu/rns/Rns_Basis.cxx', 'src/gpu/rns/Gram_Accumulator.cu', 'src/gpu/rns/Mpf_Gram.cu',
                       'src/gpu/rns/Rns_Gemm.cu', 'src/gpu/linalg/Block_Factor.cu', 'src/gpu/solver/Mpf_Ops.cu',
                       'src/gpu/solver/Gpu_Solver.cu', 'src/gpu/solver/Gpu_Collective.cu',
                       'src/gpu/linalg/Distributed_Q.cu']
        gencode = ' '.join('-gencode arch=compute_%s,code=sm_%s' % (a, a) for a in bld.env.GPU_ARCH)
        inc = '-I' + bld.env.GPU_GMP_INCLUDE if bld.env.GPU_GMP_INCLUDE else ''
        if bld.env.GPU_NCCL_DIR:
            inc += ' -DSDPB_GPU_NCCL -I' + os.path.join(bld.env.GPU_NCCL_DIR, 'include')
        gpu_headers = bld.path.ant_glob('src/gpu/**/*.hxx src/gpu/**/*.cuh')
        objs = []
        for src in gpu_sources:
            lang = '-x cu' if src.endswith('.cxx') else ''
            tgt = src + '.gpu.o'
            # waf scans no includes for a rule: every header counts, or a changed class layout leaves stale objects
            bld(rule='%s -std=c++17 -O3 %s -Xcompiler -Wall -Xcompiler -fopenmp -Xcompiler -fPIC %s %s -c ${SRC} -o ${TGT}'
                     % (bld.env.NVCC, gencode, inc, lang),
                source=src, target=tgt, deps=gpu_headers)
            objs.append(tgt)
        bld(rule='ar rcs ${TGT} ${SRC}', source=objs, target='libsdpb_gpu.a')
        # the device library is complete before anything links to it
        bld.add_group()
        lib64 = os.path.join(bld.env.CUDA_DIR, 'lib64')
        # LDFLAGS come last on the link line, after the static libraries that use them
        gpu_link = dict(ldflags=[os.path.join(bld.bldnode.abspath(), 'libsdpb_gpu.a'),
                                 os.path.join(lib64, 'libcudart_static.a'), '-L' + lib64, '-lcublasLt',
                                 '-lopenblas', '-lgomp', '-lrt', '-lpthread', '-ldl'],
                        rpath=[lib64], gpu_library=True)
    # All binaries (except for sdpb_util itself) depend also on sdpb_util
    use_packages = external_packages + ['sdpb_util']
    default_includes = ['src', 'external', bld.bldnode]

    bld.stlib(source=['src/sdpb_util/Archive_Reader.cxx',
                      'src/sdpb_util/Boost_Float.cxx',
                      'src/sdpb_util/copy_matrix.cxx',
                      'src/sdpb_util/Environment.cxx',
                      'src/sdpb_util/memory_estimates.cxx',
                      'src/sdpb_util/Memory_Tracker.cxx',
                      'src/sdpb_util/Mesh.cxx',
                      'src/sdpb_util/Proc_Meminfo.cxx',
                      'src/sdpb_util/bigint_shared_memory/blas_jobs/Blas_Job_Cost.cxx',
                      'src/sdpb_util/bigint_shared_memory/fmpz/Fmpz_BigInt.cxx',
                      'src/sdpb_util/bigint_shared_memory/fmpz/Fmpz_Matrix.cxx',
                      'src/sdpb_util/bigint_shared_memory/fmpz/Fmpz_Comb.cxx',
                      'src/sdpb_util/bigint_shared_memory/fmpz/residues.cxx',
                      'src/sdpb_util/bigint_shared_memory/fmpz/blas_residues.cxx',
                      'src/sdpb_util/block_mapping/Abstract_Block_Info.cxx',
                      'src/sdpb_util/block_mapping/Block_Mapping.cxx',
                      'src/sdpb_util/Timers/Scoped_Timer.cxx',
                      'src/sdpb_util/Timers/Timer.cxx',
                      'src/sdpb_util/Timers/Timers.cxx'],
              target='sdpb_util',
              cxxflags=default_flags,
              defines=default_defines,
              includes=default_includes,
              use=external_packages)

    sdp_solve_sources = ['src/sdp_solve/Solver_Parameters/Memory_Limit.cxx',
                         'src/sdp_solve/Solver_Parameters/Solver_Parameters.cxx',
                         'src/sdp_solve/Block_Info/Block_Info.cxx',
                         'src/sdp_solve/Block_Matrix/Block_Diagonal_Matrix.cxx',
                         'src/sdp_solve/SDP/SDP.cxx',
                         'src/sdp_solve/SDP/read_normalization.cxx',
                         'src/sdp_solve/SDP/read_objectives.cxx',
                         'src/sdp_solve/SDP/assign_bilinear_bases_dist.cxx',
                         'src/sdp_solve/SDP/read_block_data/read_block_data.cxx',
                         'src/sdp_solve/SDP/read_block_data/SDP_Block_Data.cxx',
                         'src/sdp_solve/SDP/set_bases_blocks.cxx',
                         'src/sdp_solve/SDP_Solver/save_checkpoint.cxx',
                         'src/sdp_solve/SDP_Solver/load_checkpoint/load_checkpoint.cxx',
                         'src/sdp_solve/SDP_Solver/SDP_Solver.cxx',
                         'src/sdp_solve/SDP_Solver/run/run.cxx',
                         *(['src/sdp_solve/SDP_Solver/run/gpu/Gpu_Compute.cxx'] if bld.env.SDPB_GPU else []),
                         'src/sdp_solve/SDP_Solver/run/constraint_matrix_weighted_sum.cxx',
                         'src/sdp_solve/SDP_Solver/run/compute_dual_residues_and_error.cxx',
                         'src/sdp_solve/SDP_Solver/run/compute_primal_residues_and_error_P_Ax_X.cxx',
                         'src/sdp_solve/SDP_Solver/run/compute_primal_residues_and_error_p_b_Bx.cxx',
                         'src/sdp_solve/SDP_Solver/run/compute_objectives.cxx',
                         'src/sdp_solve/SDP_Solver/run/compute_bilinear_pairings/compute_bilinear_pairings.cxx',
                         'src/sdp_solve/SDP_Solver/run/compute_bilinear_pairings/compute_A_X_inv.cxx',
                         'src/sdp_solve/SDP_Solver/run/compute_bilinear_pairings/compute_A_Y.cxx',
                         'src/sdp_solve/SDP_Solver/run/compute_feasible_and_termination.cxx',
                         'src/sdp_solve/SDP_Solver/run/print_header.cxx',
                         'src/sdp_solve/SDP_Solver/run/print_iteration.cxx',
                         'src/sdp_solve/SDP_Solver/run/step/step.cxx',
                         'src/sdp_solve/SDP_Solver/run/step/initialize_schur_complement_solver/initialize_schur_complement_solver.cxx',
                         'src/sdp_solve/SDP_Solver/run/step/initialize_schur_complement_solver/compute_schur_complement.cxx',
                         'src/sdp_solve/SDP_Solver/run/step/initialize_schur_complement_solver/compute_Q.cxx',
                         'src/sdp_solve/SDP_Solver/run/bigint_syrk/BigInt_Shared_Memory_Syrk_Context/BigInt_Shared_Memory_Syrk_Context.cxx',
                         'src/sdp_solve/SDP_Solver/run/bigint_syrk/BigInt_Shared_Memory_Syrk_Context/bigint_syrk_blas.cxx',
                         'src/sdp_solve/SDP_Solver/run/bigint_syrk/BigInt_Shared_Memory_Syrk_Context/clear_residues.cxx',
                         'src/sdp_solve/SDP_Solver/run/bigint_syrk/BigInt_Shared_Memory_Syrk_Context/compute_block_residues.cxx',
                         'src/sdp_solve/SDP_Solver/run/bigint_syrk/BigInt_Shared_Memory_Syrk_Context/get_blas_job_schedule.cxx',
                         'src/sdp_solve/SDP_Solver/run/bigint_syrk/BigInt_Shared_Memory_Syrk_Context/restore_and_reduce.cxx',
                         'src/sdp_solve/SDP_Solver/run/bigint_syrk/blas_jobs/Blas_Job.cxx',
                         'src/sdp_solve/SDP_Solver/run/bigint_syrk/blas_jobs/create_blas_job_schedule.cxx',
                         'src/sdp_solve/SDP_Solver/run/bigint_syrk/Bigint_Syrk_Config.cxx',
                         'src/sdp_solve/SDP_Solver/run/multimodular/multimodular_gemm.cxx',
                         'src/sdp_solve/SDP_Solver/run/multimodular/blocked_cholesky_trsm.cxx',
                         'src/sdp_solve/SDP_Solver/run/multimodular/Multimodular_Block_Factor.cxx',
                         'src/sdp_solve/SDP_Solver/run/step/compute_search_direction/compute_search_direction.cxx',
                         'src/sdp_solve/SDP_Solver/run/step/compute_search_direction/compute_schur_RHS.cxx',
                         'src/sdp_solve/SDP_Solver/run/step/compute_search_direction/solve_schur_complement_equation.cxx',
                         'src/sdp_solve/SDP_Solver_Terminate_Reason.cxx',
                         'src/sdp_solve/Write_Solution.cxx',
                         'src/sdp_solve/Block_Matrix/Sparse_Block_Matrix.cxx']

    bld.stlib(source=sdp_solve_sources,
              target='sdp_solve',
              cxxflags=default_flags,
              defines=default_defines,
              includes=default_includes,
              use=use_packages + ['pmp2sdp_lib'])

    sdpa_solve_sources = ['src/sdpa_solve/SDP.cxx',
                          'src/sdpa_solve/Block_Info.cxx',
                          'src/sdpa_solve/parse_sdpa/parse_sdpa.cxx',
                          'src/sdpa_solve/SDP_Solver/save_checkpoint.cxx',
                          'src/sdpa_solve/SDP_Solver/run/print_header.cxx',
                          'src/sdpa_solve/SDP_Solver/run/compute_dual_residues_and_error.cxx',
                          'src/sdpa_solve/SDP_Solver/run/step/compute_S.cxx',
                          'src/sdpa_solve/SDP_Solver/run/step/Solver_Run_Config.cxx',
                          'src/sdpa_solve/SDP_Solver/run/step/compute_search_direction.cxx',
                          'src/sdpa_solve/SDP_Solver/run/step/Initialize_P_Config.cxx',
                          'src/sdpa_solve/SDP_Solver/run/step/Initialize_P_Context.cxx',
                          'src/sdpa_solve/SDP_Solver/run/step/step.cxx',
                          'src/sdpa_solve/SDP_Solver/run/step/bigint_trmm/blas_jobs/Blas_Job.cxx',
                          'src/sdpa_solve/SDP_Solver/run/step/bigint_trmm/blas_jobs/create_blas_job_schedule.cxx',
                          'src/sdpa_solve/SDP_Solver/run/run.cxx',
                          'src/sdpa_solve/SDP_Solver/run/compute_primal_residues_and_error.cxx',
                          'src/sdpa_solve/SDP_Solver/run/print_iteration.cxx',
                          'src/sdpa_solve/SDP_Solver/run/compute_objectives.cxx',
                          'src/sdpa_solve/SDP_Solver/load_checkpoint/load_checkpoint.cxx',
                          'src/sdpa_solve/SDP_Solver/SDP_Solver.cxx']

    bld.stlib(source=sdpa_solve_sources,
              target='sdpa_solve',
              cxxflags=default_flags,
              defines=default_defines,
              includes=default_includes,
              use=use_packages + ['sdp_solve'])

    # SDPB executable
    bld.program(source=['src/sdpb/main.cxx',
                        'src/sdpb/SDPB_Parameters.cxx',
                        'src/sdpb/resident/Control_Channel.cxx',
                        'src/sdpb/resident/Resident.cxx'],
                target='sdpb',
                cxxflags=default_flags,
                defines=default_defines,
                includes=default_includes,
                use=use_packages + ['sdp_solve', 'sdpa_solve', 'mpsolve',
                                    'pmp_read', 'pmp', 'pmp2sdp_lib'],
                **gpu_link
                )

    pmp2sdp_sources = ['src/pmp2sdp/Dual_Constraint_Group/Dual_Constraint_Group.cxx',
                       'src/pmp2sdp/Dual_Constraint_Group/sample_bilinear_basis.cxx',
                       'src/pmp2sdp/Output_SDP/Output_SDP.cxx',
                       'src/pmp2sdp/write_block_data.cxx',
                       'src/pmp2sdp/write_block_info_json.cxx',
                       'src/pmp2sdp/write_objectives_json.cxx',
                       'src/pmp2sdp/write_normalization_json.cxx',
                       'src/pmp2sdp/write_sdp.cxx',
                       'src/pmp2sdp/write_control_json.cxx',
                       'src/pmp2sdp/Reuse_SDP.cxx',
                       'src/pmp2sdp/Archive_Writer.cxx',
                       'src/pmp2sdp/Archive_Entry.cxx'
                       ]

    bld.stlib(source=pmp2sdp_sources,
              target='pmp2sdp_lib',
              cxxflags=default_flags,
              defines=default_defines,
              includes=default_includes,
              use=use_packages + ['pmp'])
    bld.program(source=['src/pvm2sdp/main.cxx',
                        'src/pvm2sdp/parse_command_line.cxx'],
                target='pvm2sdp',
                cxxflags=default_flags,
                defines=default_defines,
                includes=default_includes,
                use=use_packages + ['pmp_read']
                )

    pmp_sources = ['src/pmp/Polynomial_Vector_Matrix.cxx',
                   'src/pmp/convert/sample_points.cxx',
                   'src/pmp/convert/sample_scalings.cxx',
                   'src/pmp/convert/bilinear_basis/bilinear_basis.cxx']

    bld.stlib(source=pmp_sources,
              target='pmp',
              cxxflags=default_flags,
              defines=default_defines,
              includes=default_includes,
              use=use_packages)

    pmp_read_sources = ['src/pmp_read/collect_files_expanding_nsv.cxx',
                        'src/pmp_read/PMP_File_Parse_Result.cxx',
                        'src/pmp_read/read_nsv_file_list.cxx',
                        'src/pmp_read/read_polynomial_matrix_program.cxx',
                        'src/pmp_read/read_json/read_json.cxx',
                        'src/pmp_read/read_json/Json_PMP_Parser.cxx',
                        'src/pmp_read/read_mathematica/read_mathematica.cxx',
                        'src/pmp_read/read_mathematica/parse_SDP/parse_SDP.cxx',
                        'src/pmp_read/read_mathematica/parse_SDP/parse_matrices.cxx',
                        'src/pmp_read/read_mathematica/parse_SDP/parse_number.cxx',
                        'src/pmp_read/read_mathematica/parse_SDP/parse_polynomial.cxx',
                        'src/pmp_read/read_mathematica/parse_SDP/parse_matrix/parse_matrix.cxx',
                        'src/pmp_read/read_mathematica/parse_SDP/parse_matrix/parse_damped_rational.cxx',
                        'src/pmp_read/read_xml/read_xml.cxx',
                        'src/pmp_read/read_xml/Xml_Parser/on_start_element.cxx',
                        'src/pmp_read/read_xml/Xml_Parser/on_end_element.cxx',
                        'src/pmp_read/read_xml/Xml_Parser/on_characters.cxx',
                        ]

    bld.stlib(source=pmp_read_sources,
              target='pmp_read',
              cxxflags=default_flags,
              defines=default_defines,
              includes=default_includes,
              use=use_packages + ['pmp', 'pmp2sdp_lib'])

    bld.program(source=['src/sdp2input/main.cxx'],
                target='sdp2input',
                cxxflags=default_flags,
                defines=default_defines,
                includes=default_includes,
                use=use_packages + ['pmp', 'pmp_read', 'pmp2sdp_lib']
                )

    bld.program(source=['src/pmp2sdp/main.cxx',
                        'src/pmp2sdp/Pmp2sdp_Parameters/Pmp2sdp_Parameters.cxx'
                        ],
                target='pmp2sdp',
                cxxflags=default_flags,
                defines=default_defines,
                includes=default_includes,
                use=use_packages + ['pmp', 'pmp_read', 'pmp2sdp_lib']
                )

    bld.program(source=['src/outer_limits/main.cxx',
                        'src/outer_limits/Function/eval/eval.cxx',
                        'src/outer_limits/compute_optimal/compute_optimal.cxx',
                        'src/outer_limits/compute_optimal/compute_y_transform.cxx',
                        'src/outer_limits/compute_optimal/setup_constraints.cxx',
                        'src/outer_limits/compute_optimal/find_new_points/find_new_points.cxx',
                        'src/outer_limits/compute_optimal/find_new_points/eval_summed.cxx',
                        'src/outer_limits/compute_optimal/find_new_points/get_new_points.cxx',
                        'src/outer_limits/compute_optimal/load_checkpoint/load_checkpoint.cxx',
                        'src/outer_limits/compute_optimal/load_checkpoint/Checkpoint_Parser/EndArray.cxx',
                        'src/outer_limits/compute_optimal/load_checkpoint/Checkpoint_Parser/EndObject.cxx',
                        'src/outer_limits/compute_optimal/load_checkpoint/Checkpoint_Parser/Key.cxx',
                        'src/outer_limits/compute_optimal/load_checkpoint/Checkpoint_Parser/StartArray.cxx',
                        'src/outer_limits/compute_optimal/load_checkpoint/Checkpoint_Parser/StartObject.cxx',
                        'src/outer_limits/compute_optimal/load_checkpoint/Checkpoint_Parser/String.cxx',
                        'src/outer_limits/compute_optimal/save_checkpoint.cxx',
                        'src/outer_limits/read_points/read_points.cxx',
                        'src/outer_limits/read_function_blocks/read_function_blocks.cxx',
                        'src/outer_limits/Outer_Parameters.cxx'
                        ],
                target='outer_limits',
                cxxflags=default_flags,
                defines=default_defines,
                includes=default_includes,
                use=use_packages + ['pmp_read', 'sdp_solve', 'mesh'],
                **gpu_link
                )

    bld.program(source=['src/approx_objective/main.cxx',
                        'src/approx_objective/Approx_Parameters.cxx',
                        'src/approx_objective/Axpy.cxx',
                        'src/approx_objective/setup_solver.cxx',
                        'src/approx_objective/write_solver_state.cxx',
                        'src/approx_objective/Approx_Objective/Approx_Objective.cxx',
                        'src/approx_objective/Approx_Objective/compute_dx_dy.cxx',
                        'src/approx_objective/linear_approximate_objectives.cxx',
                        'src/approx_objective/quadratic_approximate_objectives.cxx'
                        ],
                target='approx_objective',
                cxxflags=default_flags,
                defines=default_defines,
                includes=default_includes,
                use=use_packages + ['pmp_read', 'sdp_solve'],
                **gpu_link
                )

    bld.program(source=['src/pmp2functions/main.cxx',
                        'src/pmp2functions/Pmp2functions_Parameters.cxx',
                        'src/pmp2functions/write_functions.cxx'],
                target='pmp2functions',
                cxxflags=default_flags,
                defines=default_defines,
                includes=default_includes,
                use=use_packages + ['pmp_read']
                )

    bld.program(source=['src/spectrum/main.cxx',
                        'src/spectrum/handle_arguments.cxx',
                        'src/spectrum/read_c_minus_By.cxx',
                        'src/spectrum/read_pmp_info.cxx',
                        'src/spectrum/read_x.cxx',
                        'src/spectrum/set_default_parameters.cxx',
                        'src/spectrum/write_profiling.cxx',
                        'src/spectrum/compute_spectrum/compute_spectrum.cxx',
                        'src/spectrum/compute_spectrum/find_zeros.cxx',
                        'src/spectrum/compute_spectrum/mpsolve.cxx',
                        'src/spectrum/write_spectrum/write_spectrum.cxx',
                        'src/spectrum/write_spectrum/write_file.cxx'],
                target='spectrum',
                cxxflags=default_flags,
                defines=default_defines,
                includes=default_includes,
                use=use_packages + ['mpsolve']
                )

    bld.program(source=['external/catch2/catch_amalgamated.cpp',
                        'test/src/integration_tests/main.cxx',
                        'test/src/integration_tests/util/Float.cxx',
                        'test/src/integration_tests/util/diff_outer_limits.cxx',
                        'test/src/integration_tests/util/diff_sdp.cxx',
                        'test/src/integration_tests/util/diff_sdpb_out.cxx',
                        'test/src/integration_tests/util/diff_spectrum.cxx',
                        'test/src/integration_tests/util/Parse_SDP.cxx',
                        'test/src/integration_tests/util/process.cxx',
                        'test/src/integration_tests/util/Test_Case_Runner.cxx',
                        'test/src/integration_tests/cases/end-to-end.test.cxx',
                        'test/src/integration_tests/cases/outer_limits.test.cxx',
                        'test/src/integration_tests/cases/pmp2sdp.test.cxx',
                        'test/src/integration_tests/cases/SDPLIB.test.cxx',
                        'test/src/integration_tests/cases/resident.test.cxx',
                        'test/src/integration_tests/cases/sdpb.test.cxx'],
                target='integration_tests',
                install_path=None,
                cxxflags=default_flags,
                defines=default_defines + ['CATCH_AMALGAMATED_CUSTOM_MAIN'],
                use=use_packages,
                includes=default_includes + ['test/src']
                )
    bld.program(source=['external/catch2/catch_amalgamated.cpp',
                        'test/src/unit_tests/main.cxx',
                        'test/src/unit_tests/cases/LPT_scheduling.test.cxx',
                        'test/src/unit_tests/cases/Matrix_Normalizer.test.cxx',
                        'test/src/unit_tests/cases/Memory_Tracker.test.cxx',
                        'test/src/unit_tests/cases/block_data_serialization.test.cxx',
                        'test/src/unit_tests/cases/block_mapping.test.cxx',
                        'test/src/unit_tests/cases/Boost_Float.test.cxx',
                        'test/src/unit_tests/cases/blas_residues.test.cxx',
                        'test/src/unit_tests/cases/boost_serialization.test.cxx',
                        'test/src/unit_tests/cases/create_blas_job_schedule.test.cxx',
                        'test/src/unit_tests/cases/calculate_matrix_square.test.cxx',
                        'test/src/unit_tests/cases/copy_matrix.test.cxx',
                        'test/src/unit_tests/cases/json.test.cxx',
                        'test/src/unit_tests/cases/multimodular_blas.test.cxx',
                        'test/src/unit_tests/cases/pmp_sampling.test.cxx',
                        'test/src/unit_tests/cases/shared_window.test.cxx',
                        'test/src/unit_tests/cases/Sparse_Block_Matrix.test.cxx'],
                target='unit_tests',
                install_path=None,
                cxxflags=default_flags,
                defines=default_defines + ['CATCH_AMALGAMATED_CUSTOM_MAIN'],
                use=use_packages + ['pmp_read', 'pmp2sdp_lib', 'sdp_solve'],
                **gpu_link,
                includes=default_includes + ['test/src']
                )
