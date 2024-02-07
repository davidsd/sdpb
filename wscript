import os, subprocess


def options(opt):
    opt.load(['compiler_cxx', 'gnu_dirs'])
    opt.load(['cxx17', 'boost', 'gmpxx', 'mpfr', 'elemental', 'libxml2', 'rapidjson', 'libarchive'],
             tooldir='./waf-tools')


def configure(conf):
    if not 'CXX' in os.environ or os.environ['CXX'] == 'g++' or os.environ['CXX'] == 'icpc':
        conf.environ['CXX'] = 'mpicxx'

    conf.load(['compiler_cxx', 'gnu_dirs', 'cxx17', 'boost', 'gmpxx', 'mpfr',
               'elemental', 'libxml2', 'rapidjson', 'libarchive'])
    conf.load('clang_compilation_database', tooldir='./waf-tools')

    conf.env.git_version = subprocess.check_output('git describe --tags --always --dirty', universal_newlines=True,
                                                   shell=True).rstrip()


def build(bld):
    default_flags = ['-Wall', '-Wextra', '-Werror=return-type', '-O3']
    default_defines = ['OMPI_SKIP_MPICXX', 'SDPB_VERSION_STRING="' + bld.env.git_version + '"']
    use_packages = ['cxx17', 'gmpxx', 'mpfr', 'boost', 'elemental', 'libxml2', 'rapidjson', 'libarchive', 'sdpb_util']
    default_includes = ['src', 'external']

    bld.stlib(source=['src/sdpb_util/copy_matrix.cxx',
                      'src/sdpb_util/Environment.cxx',
                      'src/sdpb_util/Mesh.cxx',
                      'src/sdpb_util/Proc_Meminfo.cxx',
                      'src/sdpb_util/Timers/Scoped_Timer.cxx',
                      'src/sdpb_util/Timers/Timer.cxx',
                      'src/sdpb_util/Timers/Timers.cxx'],
              target='sdpb_util',
              cxxflags=default_flags,
              defines=default_defines,
              includes=default_includes,
              use=['cxx17', 'gmpxx', 'boost', 'elemental'])

    sdp_solve_sources = ['src/sdp_solve/Solver_Parameters/Solver_Parameters.cxx',
                         'src/sdp_solve/Solver_Parameters/ostream.cxx',
                         'src/sdp_solve/Solver_Parameters/to_property_tree.cxx',
                         'src/sdp_solve/Archive_Reader/Archive_Reader.cxx',
                         'src/sdp_solve/Archive_Reader/underflow.cxx',
                         'src/sdp_solve/Block_Info/Block_Info.cxx',
                         'src/sdp_solve/Block_Info/read_block_info.cxx',
                         'src/sdp_solve/Block_Info/read_block_costs.cxx',
                         'src/sdp_solve/Block_Info/allocate_blocks.cxx',
                         'src/sdp_solve/SDP/SDP/SDP.cxx',
                         'src/sdp_solve/SDP/SDP/read_normalization.cxx',
                         'src/sdp_solve/SDP/SDP/read_objectives.cxx',
                         'src/sdp_solve/SDP/SDP/assign_bilinear_bases_dist.cxx',
                         'src/sdp_solve/SDP/SDP/read_block_data/read_block_data.cxx',
                         'src/sdp_solve/SDP/SDP/read_block_data/SDP_Block_Data.cxx',
                         'src/sdp_solve/SDP/SDP/set_bases_blocks.cxx',
                         'src/sdp_solve/SDP_Solver/save_checkpoint.cxx',
                         'src/sdp_solve/SDP_Solver/load_checkpoint/load_checkpoint.cxx',
                         'src/sdp_solve/SDP_Solver/load_checkpoint/load_binary_checkpoint.cxx',
                         'src/sdp_solve/SDP_Solver/load_checkpoint/load_text_checkpoint.cxx',
                         'src/sdp_solve/SDP_Solver/SDP_Solver.cxx',
                         'src/sdp_solve/SDP_Solver/run/run.cxx',
                         'src/sdp_solve/SDP_Solver/run/cholesky_decomposition.cxx',
                         'src/sdp_solve/SDP_Solver/run/constraint_matrix_weighted_sum.cxx',
                         'src/sdp_solve/SDP_Solver/run/compute_dual_residues_and_error.cxx',
                         'src/sdp_solve/SDP_Solver/run/compute_primal_residues_and_error_P_Ax_X.cxx',
                         'src/sdp_solve/SDP_Solver/run/compute_primal_residues_and_error_p_b_Bx.cxx',
                         'src/sdp_solve/SDP_Solver/run/compute_objectives/compute_objectives.cxx',
                         'src/sdp_solve/SDP_Solver/run/compute_objectives/dot.cxx',
                         'src/sdp_solve/SDP_Solver/run/compute_bilinear_pairings/compute_bilinear_pairings.cxx',
                         'src/sdp_solve/SDP_Solver/run/compute_bilinear_pairings/compute_A_X_inv.cxx',
                         'src/sdp_solve/SDP_Solver/run/compute_bilinear_pairings/compute_A_Y.cxx',
                         'src/sdp_solve/SDP_Solver/run/compute_feasible_and_termination.cxx',
                         'src/sdp_solve/SDP_Solver/run/print_header.cxx',
                         'src/sdp_solve/SDP_Solver/run/print_iteration.cxx',
                         'src/sdp_solve/SDP_Solver/run/step/step.cxx',
                         'src/sdp_solve/SDP_Solver/run/step/initialize_schur_complement_solver/initialize_schur_complement_solver.cxx',
                         'src/sdp_solve/SDP_Solver/run/step/initialize_schur_complement_solver/compute_schur_complement.cxx',
                         'src/sdp_solve/SDP_Solver/run/step/initialize_schur_complement_solver/initialize_Q_group.cxx',
                         'src/sdp_solve/SDP_Solver/run/step/initialize_schur_complement_solver/synchronize_Q.cxx',
                         'src/sdp_solve/SDP_Solver/run/step/compute_search_direction/compute_search_direction.cxx',
                         'src/sdp_solve/SDP_Solver/run/step/compute_search_direction/cholesky_solve.cxx',
                         'src/sdp_solve/SDP_Solver/run/step/compute_search_direction/compute_schur_RHS.cxx',
                         'src/sdp_solve/SDP_Solver/run/step/compute_search_direction/scale_multiply_add.cxx',
                         'src/sdp_solve/SDP_Solver/run/step/compute_search_direction/solve_schur_complement_equation.cxx',
                         'src/sdp_solve/SDP_Solver/run/step/predictor_centering_parameter.cxx',
                         'src/sdp_solve/SDP_Solver/run/step/corrector_centering_parameter/corrector_centering_parameter.cxx',
                         'src/sdp_solve/SDP_Solver/run/step/corrector_centering_parameter/frobenius_product_of_sums.cxx',
                         'src/sdp_solve/SDP_Solver/run/step/frobenius_product_symmetric.cxx',
                         'src/sdp_solve/SDP_Solver/run/step/step_length/step_length.cxx',
                         'src/sdp_solve/SDP_Solver/run/step/step_length/min_eigenvalue.cxx',
                         'src/sdp_solve/SDP_Solver/run/step/step_length/lower_triangular_inverse_congruence.cxx',
                         'src/sdp_solve/SDP_Solver_Terminate_Reason/ostream.cxx',
                         'src/sdp_solve/lower_triangular_transpose_solve.cxx',
                         'src/sdp_solve/Block_Diagonal_Matrix/ostream.cxx',
                         'src/sdp_solve/Write_Solution.cxx']

    bld.stlib(source=sdp_solve_sources,
              target='sdp_solve',
              cxxflags=default_flags,
              defines=default_defines,
              includes=default_includes,
              use=use_packages + ['pmp2sdp_lib'])

    # SDPB executable
    bld.program(source=['src/sdpb/main.cxx',
                        'src/sdpb/solve.cxx',
                        'src/sdpb/write_timing.cxx',
                        'src/sdpb/SDPB_Parameters/SDPB_Parameters.cxx',
                        'src/sdpb/SDPB_Parameters/to_property_tree.cxx',
                        'src/sdpb/SDPB_Parameters/ostream.cxx',
                        'src/sdpb/save_solution.cxx'],
                target='sdpb',
                cxxflags=default_flags,
                defines=default_defines,
                includes=default_includes,
                use=use_packages + ['sdp_solve']
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
                       'src/pmp2sdp/Archive_Writer/Archive_Writer.cxx',
                       'src/pmp2sdp/Archive_Writer/write_entry.cxx',
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
                   'src/pmp/convert/bilinear_basis/bilinear_basis.cxx',
                   'src/pmp/convert/bilinear_basis/precompute/precompute.cxx',
                   'src/pmp/convert/bilinear_basis/precompute/integral.cxx',
                   'src/pmp/convert/bilinear_basis/bilinear_form/bilinear_form.cxx',
                   'src/pmp/convert/bilinear_basis/bilinear_form/rest.cxx',
                   'src/pmp/convert/bilinear_basis/bilinear_form/dExp.cxx',
                   'src/pmp/convert/bilinear_basis/bilinear_form/derivative.cxx',
                   'src/pmp/convert/bilinear_basis/bilinear_form/operator_plus_set_Derivative_Term.cxx']

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
                        'src/outer_limits/power_prefactor.cxx',
                        'src/outer_limits/poles_prefactor.cxx',
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
                        'src/outer_limits/Outer_Parameters/Outer_Parameters.cxx',
                        'src/outer_limits/Outer_Parameters/to_property_tree.cxx',
                        'src/outer_limits/Outer_Parameters/ostream.cxx'
                        ],
                target='outer_limits',
                cxxflags=default_flags,
                defines=default_defines,
                includes=default_includes,
                use=use_packages + ['pmp_read', 'sdp_solve', 'mesh']
                )

    bld.program(source=['src/approx_objective/main.cxx',
                        'src/approx_objective/Approx_Parameters/Approx_Parameters.cxx',
                        'src/approx_objective/Approx_Parameters/ostream.cxx',
                        'src/approx_objective/Axpy.cxx',
                        'src/approx_objective/setup_solver.cxx',
                        'src/approx_objective/write_solver_state.cxx',
                        'src/approx_objective/Approx_Objective/Approx_Objective/Approx_Objective.cxx',
                        'src/approx_objective/Approx_Objective/Approx_Objective/compute_dx_dy.cxx',
                        'src/approx_objective/linear_approximate_objectives.cxx',
                        'src/approx_objective/quadratic_approximate_objectives.cxx'
                        ],
                target='approx_objective',
                cxxflags=default_flags,
                defines=default_defines,
                includes=default_includes,
                use=use_packages + ['pmp_read', 'sdp_solve']
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
                        'src/spectrum/read_x.cxx',
                        'src/spectrum/compute_spectrum.cxx',
                        'src/spectrum/compute_lambda.cxx',
                        'src/spectrum/eval_summed.cxx',
                        'src/spectrum/get_zeros.cxx',
                        'src/spectrum/write_spectrum/write_spectrum.cxx',
                        'src/spectrum/write_spectrum/write_file.cxx'],
                target='spectrum',
                cxxflags=default_flags,
                defines=default_defines,
                includes=default_includes,
                use=use_packages + ['pmp_read', 'sdp_solve', 'pmp2sdp_lib', 'mesh']
                )

    bld.program(source=['external/catch2/catch_amalgamated.cpp',
                        'test/src/integration_tests/main.cxx',
                        'test/src/integration_tests/util/Float.cxx',
                        'test/src/integration_tests/util/diff_outer_limits.cxx',
                        'test/src/integration_tests/util/diff_sdp.cxx',
                        'test/src/integration_tests/util/diff_sdpb_out.cxx',
                        'test/src/integration_tests/util/diff_spectrum.cxx',
                        'test/src/integration_tests/util/Test_Case_Runner.cxx',
                        'test/src/integration_tests/cases/end-to-end.test.cxx',
                        'test/src/integration_tests/cases/outer_limits.test.cxx',
                        'test/src/integration_tests/cases/pmp2sdp.test.cxx',
                        'test/src/integration_tests/cases/sdpb.test.cxx',
                        'test/src/integration_tests/cases/spectrum.test.cxx'],
                target='integration_tests',
                install_path=None,
                cxxflags=default_flags,
                defines=default_defines + ['CATCH_AMALGAMATED_CUSTOM_MAIN'],
                use=use_packages,
                includes=default_includes + ['test/src']
                )
    bld.program(source=['external/catch2/catch_amalgamated.cpp',
                        'test/src/unit_tests/main.cxx',
                        'test/src/unit_tests/cases/block_data_serialization.test.cxx',
                        'test/src/unit_tests/cases/block_mapping.test.cxx',
                        'test/src/unit_tests/cases/Boost_Float.test.cxx',
                        'test/src/unit_tests/cases/boost_serialization.test.cxx',
                        'test/src/unit_tests/cases/copy_matrix.test.cxx',
                        'test/src/unit_tests/cases/json.test.cxx',
                        'test/src/unit_tests/cases/shared_window.test.cxx'],
                target='unit_tests',
                cxxflags=default_flags,
                defines=default_defines + ['CATCH_AMALGAMATED_CUSTOM_MAIN'],
                use=use_packages + ['pmp_read', 'pmp2sdp_lib', 'sdp_solve'],
                includes=default_includes + ['test/src']
                )
