import os, subprocess

def options(opt):
    opt.load(['compiler_cxx','gnu_dirs','cxx17','boost','gmpxx','mpfr',
              'elemental','libxml2', 'rapidjson', 'libarchive'])

def configure(conf):
    if not 'CXX' in os.environ or os.environ['CXX']=='g++' or os.environ['CXX']=='icpc':
        conf.environ['CXX']='mpicxx'

    conf.load(['compiler_cxx','gnu_dirs','cxx17','boost','gmpxx','mpfr',
               'elemental','libxml2', 'rapidjson', 'libarchive'])
    conf.load('clang_compilation_database')

    conf.env.git_version=subprocess.check_output('git describe --dirty', universal_newlines=True, shell=True).rstrip()
    
def build(bld):
    # optimized building:
    # default_flags=['-Wall', '-Wextra', '-O3', '-DOMPI_SKIP_MPICXX', '-D SDPB_VERSION_STRING="' + bld.env.git_version + '"']
    
    # default_flags=['-Wall', '-Wextra', '-O3', '-g', '-DOMPI_SKIP_MPICXX', '-D SDPB_VERSION_STRING="' + bld.env.git_version + '"']
    default_flags=['-Wall', '-Wextra', '-g', '-DOMPI_SKIP_MPICXX', '-D SDPB_VERSION_STRING="' + bld.env.git_version + '"']
    use_packages=['cxx17','boost','gmpxx','mpfr','elemental','libxml2', 'rapidjson', 'libarchive']

    sdp_solve_sources=['src/sdp_solve/Solver_Parameters/Solver_Parameters.cxx',
                       'src/sdp_solve/Solver_Parameters/ostream.cxx',
                       'src/sdp_solve/Solver_Parameters/to_property_tree.cxx',
                       'src/sdp_solve/Archive_Reader/Archive_Reader.cxx',
                       'src/sdp_solve/Archive_Reader/underflow.cxx',
                       'src/sdp_solve/Block_Info/Block_Info.cxx',
                       'src/sdp_solve/Block_Info/read_block_info.cxx',
                       'src/sdp_solve/Block_Info/read_block_costs.cxx',
                       'src/sdp_solve/Block_Info/allocate_blocks/allocate_blocks.cxx',
                       'src/sdp_solve/Block_Info/allocate_blocks/compute_block_grid_mapping.cxx',
                       'src/sdp_solve/SDP/SDP/SDP.cxx',
                       'src/sdp_solve/SDP/SDP/read_objectives.cxx',
                       'src/sdp_solve/SDP/SDP/assign_bilinear_bases_dist.cxx',
                       'src/sdp_solve/SDP/SDP/read_blocks/read_blocks.cxx',
                       'src/sdp_solve/SDP/SDP/read_blocks/read_block_stream/read_block_stream.cxx',
                       'src/sdp_solve/SDP/SDP/read_blocks/read_block_stream/Block_Parser/EndArray.cxx',
                       'src/sdp_solve/SDP/SDP/read_blocks/read_block_stream/Block_Parser/Key.cxx',
                       'src/sdp_solve/SDP/SDP/read_blocks/read_block_stream/Block_Parser/StartArray.cxx',
                       'src/sdp_solve/SDP/SDP/read_blocks/read_block_stream/Block_Parser/String.cxx',
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
              use=use_packages + ['sdp_read'])

    

    dynamical_solve_sources=['src/dynamical_solve/Dynamical_Solver_Parameters/Dynamical_Solver_Parameters.cxx',
                       'src/dynamical_solve/Dynamical_Solver_Parameters/ostream.cxx',
                       'src/dynamical_solve/Dynamical_Solver_Parameters/to_property_tree.cxx',    
                       'src/dynamical_solve/Dynamical_Solver_Terminate_Reason/ostream.cxx', 
                       'src/dynamical_solve/Dynamical_Solver/save_checkpoint.cxx',
                       'src/dynamical_solve/Dynamical_Solver/load_checkpoint/load_checkpoint.cxx',
                       'src/dynamical_solve/Dynamical_Solver/load_checkpoint/load_binary_checkpoint.cxx',
                       'src/dynamical_solve/Dynamical_Solver/load_checkpoint/load_text_checkpoint.cxx',
                       'src/dynamical_solve/Dynamical_Solver/Dynamical_Solver.cxx',
                       'src/dynamical_solve/Dynamical_Solver/dynamical_navigator/mixed_hess.cxx',
                       'src/dynamical_solve/Dynamical_Solver/dynamical_navigator/external_grad_hess.cxx',
                       'src/dynamical_solve/Dynamical_Solver/dynamical_navigator/compute_lag.cxx',
                       #'src/dynamical_solve/Dynamical_Solver/dynamical_navigator/dynamical_step_robust.cxx',
                       'src/dynamical_solve/Dynamical_Solver/dynamical_navigator/BFGS.cxx',
                       'src/dynamical_solve/Dynamical_Solver/dynamical_navigator/dynamical_step_BFGS.cxx',
                       'src/dynamical_solve/Dynamical_Solver/dynamical_navigator/dynamical_step_BFGS_func.cxx',
                       'src/dynamical_solve/Dynamical_Solver/dynamical_navigator/dynamical_step_external_corrector.cxx',
                       'src/dynamical_solve/Dynamical_Solver/dynamical_navigator/internal_search_direction.cxx', 
                       'src/dynamical_solve/Dynamical_Solver/dynamical_navigator/compute_update_sdp.cxx',
                       'src/dynamical_solve/Dynamical_Solver/dynamical_navigator/compute_search_direction.cxx',
                       'src/dynamical_solve/Dynamical_Solver/run_dynamical.cxx',
                       'src/dynamical_solve/Dynamical_Solver/run/print_header_dynamical.cxx',
                       'src/dynamical_solve/Dynamical_Solver/run/print_iteration.cxx',
                       'src/dynamical_solve/Dynamical_Solver/dynamical_navigator/compute_R_error.cxx',
                       'src/dynamical_solve/Dynamical_Solver/run/compute_feasible_and_termination.cxx',
                       'src/dynamical_solve/Dynamical_Solver/dynamical_navigator/cubic_approx.cxx',
                       'src/approx_objective/Approx_Parameters/Approx_Parameters.cxx',
                       'src/approx_objective/Axpy.cxx',
                       'src/approx_objective/Approx_Objective/Approx_Objective/Approx_Objective.cxx',
                       'src/approx_objective/Approx_Objective/Approx_Objective/compute_dx_dy.cxx',
                       'src/dynamical_solve/Dynamical_Solver/dynamical_navigator/newton_trust.cxx',
                       'src/approx_objective/write_solver_state.cxx']

    bld.stlib(source=dynamical_solve_sources,
          target='dynamical_solve',
          cxxflags=default_flags,
          use=use_packages + ['sdp_solve', 'sdp_read'])

                      
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
                use=use_packages + ['sdp_solve']
                )

    sdp_convert_sources=['src/sdp_convert/Dual_Constraint_Group/Dual_Constraint_Group/Dual_Constraint_Group.cxx',
                         'src/sdp_convert/Dual_Constraint_Group/Dual_Constraint_Group/sample_bilinear_basis.cxx',
                         'src/sdp_convert/write_objectives.cxx',
                         'src/sdp_convert/write_bilinear_bases.cxx',
                         'src/sdp_convert/write_blocks.cxx',
                         'src/sdp_convert/write_primal_objective_c.cxx',
                         'src/sdp_convert/write_free_var_matrix.cxx',
                         'src/sdp_convert/write_sdpb_input_files.cxx',
                         'src/sdp_convert/write_control.cxx',
                         'src/sdp_convert/Archive_Writer/Archive_Writer.cxx',
                         'src/sdp_convert/Archive_Writer/write_entry.cxx',
                         'src/sdp_convert/Archive_Entry.cxx'
                         ]

    bld.stlib(source=sdp_convert_sources,
              target='sdp_convert',
              cxxflags=default_flags,
              use=use_packages)


    bld.program(source=['src/pvm2sdp/main.cxx',
                        'src/pvm2sdp/parse_command_line.cxx',
                        'src/pvm2sdp/read_input_files/read_input_files.cxx',
                        'src/pvm2sdp/read_input_files/read_xml_input/read_xml_input.cxx',
                        'src/pvm2sdp/read_input_files/read_xml_input/Input_Parser/on_start_element.cxx',
                        'src/pvm2sdp/read_input_files/read_xml_input/Input_Parser/on_end_element.cxx',
                        'src/pvm2sdp/read_input_files/read_xml_input/Input_Parser/on_characters.cxx'],
                target='pvm2sdp',
                cxxflags=default_flags,
                use=use_packages + ['sdp_read']
                )

    sdp_read_sources=['src/sdp_read/read_input/read_input.cxx',
                      'src/sdp_read/read_input/read_json/read_json.cxx',
                      'src/sdp_read/read_input/read_json/Positive_Matrix_With_Prefactor_State/json_key.cxx',
                      'src/sdp_read/read_input/read_json/Positive_Matrix_With_Prefactor_State/json_string.cxx',
                      'src/sdp_read/read_input/read_json/Positive_Matrix_With_Prefactor_State/json_start_array.cxx',
                      'src/sdp_read/read_input/read_json/Positive_Matrix_With_Prefactor_State/json_end_array.cxx',
                      'src/sdp_read/read_input/read_json/Positive_Matrix_With_Prefactor_State/json_start_object.cxx',
                      'src/sdp_read/read_input/read_json/Positive_Matrix_With_Prefactor_State/json_end_object.cxx',
                      'src/sdp_read/read_input/read_json/Damped_Rational_State/json_key.cxx',
                      'src/sdp_read/read_input/read_json/Damped_Rational_State/json_string.cxx',
                      'src/sdp_read/read_input/read_json/Damped_Rational_State/json_start_array.cxx',
                      'src/sdp_read/read_input/read_json/Damped_Rational_State/json_end_array.cxx',
                      'src/sdp_read/read_input/read_json/Damped_Rational_State/json_start_object.cxx',
                      'src/sdp_read/read_input/read_json/Damped_Rational_State/json_end_object.cxx',
                      'src/sdp_read/read_input/read_json/JSON_Parser/Key.cxx',
                      'src/sdp_read/read_input/read_json/JSON_Parser/String.cxx',
                      'src/sdp_read/read_input/read_json/JSON_Parser/StartArray.cxx',
                      'src/sdp_read/read_input/read_json/JSON_Parser/EndArray.cxx',
                      'src/sdp_read/read_input/read_json/JSON_Parser/StartObject.cxx',
                      'src/sdp_read/read_input/read_json/JSON_Parser/EndObject.cxx',
                      'src/sdp_read/read_input/read_mathematica/read_mathematica.cxx',
                      'src/sdp_read/read_input/read_mathematica/parse_SDP/parse_SDP.cxx',
                      'src/sdp_read/read_input/read_mathematica/parse_SDP/parse_matrices.cxx',
                      'src/sdp_read/read_input/read_mathematica/parse_SDP/parse_number.cxx',
                      'src/sdp_read/read_input/read_mathematica/parse_SDP/parse_polynomial.cxx',
                      'src/sdp_read/read_input/read_mathematica/parse_SDP/parse_matrix/parse_matrix.cxx',
                      'src/sdp_read/read_input/read_mathematica/parse_SDP/parse_matrix/parse_damped_rational.cxx',
                      'src/sdp_read/read_pvm_input/read_pvm_input.cxx',
                      'src/sdp_read/read_pvm_input/read_xml_input/read_xml_input.cxx',
                      'src/sdp_read/read_pvm_input/read_xml_input/Input_Parser/on_start_element.cxx',
                      'src/sdp_read/read_pvm_input/read_xml_input/Input_Parser/on_end_element.cxx',
                      'src/sdp_read/read_pvm_input/read_xml_input/Input_Parser/on_characters.cxx',
                      'src/sdp_read/read_file_list.cxx',
                      'src/sdp_read/sample_points.cxx',
                      'src/sdp_read/sample_scalings.cxx'
                      ]

    bld.stlib(source=sdp_read_sources,
              target='sdp_read',
              cxxflags=default_flags,
              use=use_packages + ['sdp_convert'])

    bld.program(source=['src/sdp2input/main.cxx',
                        'src/sdp2input/write_output/write_output.cxx',
                        'src/sdp2input/write_output/bilinear_basis/bilinear_basis.cxx',
                        'src/sdp2input/write_output/bilinear_basis/precompute/precompute.cxx',
                        'src/sdp2input/write_output/bilinear_basis/precompute/integral.cxx',
                        'src/sdp2input/write_output/bilinear_basis/bilinear_form/bilinear_form.cxx',
                        'src/sdp2input/write_output/bilinear_basis/bilinear_form/rest.cxx',
                        'src/sdp2input/write_output/bilinear_basis/bilinear_form/dExp.cxx',
                        'src/sdp2input/write_output/bilinear_basis/bilinear_form/derivative.cxx',
                        'src/sdp2input/write_output/bilinear_basis/bilinear_form/operator_plus_set_Derivative_Term.cxx'],
                target='sdp2input',
                cxxflags=default_flags,
                use=use_packages + ['sdp_read']
                )

    bld.stlib(source=['src/Mesh/Mesh.cxx', 'src/Mesh/ostream.cxx'],
              target='mesh',
              cxxflags=default_flags,
              use=use_packages)

    bld.program(source=['src/outer_limits/main.cxx',
                        'src/outer_limits/power_prefactor.cxx',
                        'src/outer_limits/poles_prefactor.cxx',
                        'src/outer_limits/Function/eval/eval.cxx',
                        'src/outer_limits/compute_optimal/compute_optimal.cxx',
                        'src/outer_limits/compute_optimal/compute_y_transform.cxx',
                        'src/outer_limits/compute_optimal/setup_constraints.cxx',
                        'src/outer_limits/compute_optimal/copy_matrix.cxx',
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
                        'src/outer_limits/read_points/read_points_json/read_points_json.cxx',
                        'src/outer_limits/read_points/read_points_json/Points_Parser/EndArray.cxx',
                        'src/outer_limits/read_points/read_points_json/Points_Parser/EndObject.cxx',
                        'src/outer_limits/read_points/read_points_json/Points_Parser/Key.cxx',
                        'src/outer_limits/read_points/read_points_json/Points_Parser/StartArray.cxx',
                        'src/outer_limits/read_points/read_points_json/Points_Parser/StartObject.cxx',
                        'src/outer_limits/read_points/read_points_json/Points_Parser/String.cxx',
                        'src/outer_limits/read_function_blocks/read_function_blocks.cxx',
                        'src/outer_limits/read_function_blocks/read_json/Function_State/json_end_array.cxx',
                        'src/outer_limits/read_function_blocks/read_json/Function_State/json_end_object.cxx',
                        'src/outer_limits/read_function_blocks/read_json/Function_State/json_key.cxx',
                        'src/outer_limits/read_function_blocks/read_json/Function_State/json_start_array.cxx',
                        'src/outer_limits/read_function_blocks/read_json/Function_State/json_start_object.cxx',
                        'src/outer_limits/read_function_blocks/read_json/Function_State/json_string.cxx',
                        'src/outer_limits/read_function_blocks/read_json/Function_Blocks_Parser/EndArray.cxx',
                        'src/outer_limits/read_function_blocks/read_json/Function_Blocks_Parser/EndObject.cxx',
                        'src/outer_limits/read_function_blocks/read_json/Function_Blocks_Parser/Key.cxx',
                        'src/outer_limits/read_function_blocks/read_json/Function_Blocks_Parser/StartArray.cxx',
                        'src/outer_limits/read_function_blocks/read_json/Function_Blocks_Parser/StartObject.cxx',
                        'src/outer_limits/read_function_blocks/read_json/Function_Blocks_Parser/String.cxx',
                        'src/outer_limits/read_function_blocks/read_json/read_json.cxx',
                        'src/outer_limits/Outer_Parameters/Outer_Parameters.cxx',
                        'src/outer_limits/Outer_Parameters/to_property_tree.cxx',
                        'src/outer_limits/Outer_Parameters/ostream.cxx'
                        ],
                target='outer_limits',
                cxxflags=default_flags,
                use=use_packages + ['sdp_read','sdp_solve', 'mesh']
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
                use=use_packages + ['sdp_read','sdp_solve']
                )
    
    bld.program(source=['src/pvm2functions/main.cxx',
                        'src/pvm2functions/parse_command_line.cxx',
                        'src/pvm2functions/write_functions.cxx'],
                target='pvm2functions',
                cxxflags=default_flags,
                use=use_packages + ['sdp_read']
                )

    bld.program(source=['src/sdp2functions/main.cxx',
                        'src/sdp2functions/write_functions.cxx'],
                target='sdp2functions',
                cxxflags=default_flags,
                use=use_packages + ['sdp_read']
                )

    bld.program(source=['src/spectrum/main.cxx',
                        'src/spectrum/handle_arguments.cxx',
                        'src/spectrum/read_x.cxx',
                        'src/spectrum/compute_spectrum_pmp.cxx',
                        'src/spectrum/compute_spectrum_pvm.cxx',
                        'src/spectrum/compute_lambda.cxx',
                        'src/spectrum/eval_summed.cxx',
                        'src/spectrum/get_zeros.cxx',
                        'src/spectrum/write_spectrum/write_spectrum.cxx',
                        'src/spectrum/write_spectrum/write_file.cxx'],
                target='spectrum',
                cxxflags=default_flags,
                use=use_packages + ['sdp_read', 'sdp_solve', 'sdp_convert', 'mesh']
                )
    

    # Dynamically Navigated SDP executable
    bld.program(source=['src/dynamical_sdp/main.cxx',
                        'src/dynamical_sdp/solve.cxx',
                        'src/dynamical_sdp/write_timing.cxx',
                        'src/dynamical_sdp/Dynamical_Parameters/Dynamical_Parameters.cxx',
                        'src/dynamical_sdp/Dynamical_Parameters/to_property_tree.cxx',
                        'src/dynamical_sdp/Dynamical_Parameters/ostream.cxx',
                        'src/dynamical_sdp/save_solution.cxx'],
                target='dynamical_sdp_V15G',
                cxxflags=default_flags,
                use=use_packages + ['sdp_read','sdp_solve', 'dynamical_solve']
                )


    # Test NewtonRegion
    #bld.program(source=['src/dynamical_solve/Dynamical_Solver/dynamical_navigator/newton_trust.cxx'],
    #            target='test_newton',
    #            cxxflags=default_flags,
    #            use=use_packages + ['sdp_solve', 'dynamical_solve']
    #            )
