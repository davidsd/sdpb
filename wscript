import os

def options(opt):
    opt.load(['compiler_cxx','gnu_dirs','cxx14','boost','gmpxx','mpfr',
              'elemental','libxml2'])

def configure(conf):
    if not 'CXX' in os.environ or os.environ['CXX']=='g++' or os.environ['CXX']=='icpc':
        conf.environ['CXX']='mpicxx'

    conf.load(['compiler_cxx','gnu_dirs','cxx14','boost','gmpxx','mpfr',
               'elemental','libxml2'])

def build(bld):
    default_flags=['-Wall', '-Wextra', '-O3']
    # default_flags=['-Wall', '-Wextra', '-g']
    use_packages=['cxx14','boost','gmpxx','mpfr','elemental','libxml2']
    
    # Main executable
    bld.program(source=['src/sdpb/main.cxx',
                        'src/sdpb/SDP_Solver_Parameters/ostream.cxx',
                        'src/sdpb/solve/solve.cxx',
                        'src/compute_block_grid_mapping.cxx',
                        'src/sdpb/Block_Info/Block_Info.cxx',
                        'src/sdpb/write_timing.cxx',
                        'src/sdpb/solve/SDP/SDP/SDP.cxx',
                        'src/sdpb/solve/SDP/SDP/read_objectives.cxx',
                        'src/sdpb/solve/SDP/SDP/read_bilinear_bases.cxx',
                        'src/sdpb/solve/SDP/SDP/read_primal_objective_c.cxx',
                        'src/sdpb/solve/SDP/SDP/read_free_var_matrix.cxx',
                        'src/sdpb/solve/SDP_Solver/save_solution.cxx',
                        'src/sdpb/solve/SDP_Solver/save_checkpoint.cxx',
                        'src/sdpb/solve/SDP_Solver/load_checkpoint.cxx',
                        'src/sdpb/solve/SDP_Solver/SDP_Solver.cxx',
                        'src/sdpb/solve/SDP_Solver/run/run.cxx',
                        'src/sdpb/solve/SDP_Solver/run/constraint_matrix_weighted_sum.cxx',
                        'src/sdpb/solve/SDP_Solver/run/compute_dual_residues_and_error.cxx',
                        'src/sdpb/solve/SDP_Solver/run/compute_primal_residues_and_error.cxx',
                        'src/sdpb/solve/SDP_Solver/run/compute_objectives/compute_objectives.cxx',
                        'src/sdpb/solve/SDP_Solver/run/compute_objectives/dot.cxx',
                        'src/sdpb/solve/SDP_Solver/run/compute_bilinear_pairings/compute_bilinear_pairings.cxx',
                        'src/sdpb/solve/SDP_Solver/run/compute_bilinear_pairings/compute_bilinear_pairings_X_inv.cxx',
                        'src/sdpb/solve/SDP_Solver/run/compute_bilinear_pairings/compute_bilinear_pairings_Y.cxx',
                        'src/sdpb/solve/SDP_Solver/run/compute_feasible_and_termination.cxx',
                        'src/sdpb/solve/SDP_Solver/run/print_header.cxx',
                        'src/sdpb/solve/SDP_Solver/run/step/step.cxx',
                        'src/sdpb/solve/SDP_Solver/run/step/print_iteration.cxx',
                        'src/sdpb/solve/SDP_Solver/run/step/initialize_schur_complement_solver/initialize_schur_complement_solver.cxx',
                        'src/sdpb/solve/SDP_Solver/run/step/initialize_schur_complement_solver/compute_schur_complement.cxx',
                        'src/sdpb/solve/SDP_Solver/run/step/initialize_schur_complement_solver/synchronize/synchronize.cxx',
                        'src/sdpb/solve/SDP_Solver/run/step/initialize_schur_complement_solver/synchronize/reduce_and_scatter.cxx',
                        'src/sdpb/solve/SDP_Solver/run/step/compute_search_direction/compute_search_direction.cxx',
                        'src/sdpb/solve/SDP_Solver/run/step/compute_search_direction/compute_schur_RHS.cxx',
                        'src/sdpb/solve/SDP_Solver/run/step/compute_search_direction/solve_schur_complement_equation.cxx',
                        'src/sdpb/solve/SDP_Solver/run/step/predictor_centering_parameter.cxx',
                        'src/sdpb/solve/SDP_Solver/run/step/corrector_centering_parameter.cxx',
                        'src/sdpb/solve/SDP_Solver/run/step/step_length.cxx',
                        'src/sdpb/solve/SDP_Solver_Terminate_Reason/ostream.cxx',
                        'src/sdpb/solve/lower_triangular_transpose_solve.cxx',
                        'src/sdpb/solve/Block_Diagonal_Matrix/scale_multiply_add.cxx',
                        'src/sdpb/solve/Block_Diagonal_Matrix/cholesky_solve.cxx',
                        'src/sdpb/solve/Block_Diagonal_Matrix/cholesky_decomposition.cxx',
                        'src/sdpb/solve/Block_Diagonal_Matrix/frobenius_product_of_sums.cxx',
                        'src/sdpb/solve/Block_Diagonal_Matrix/frobenius_product_symmetric.cxx',
                        'src/sdpb/solve/Block_Diagonal_Matrix/lower_triangular_inverse_congruence.cxx',
                        'src/sdpb/solve/Block_Diagonal_Matrix/min_eigenvalue.cxx',
                        'src/sdpb/solve/Block_Diagonal_Matrix/ostream.cxx'],
                target='sdpb',
                cxxflags=default_flags,
                use=use_packages
                )
    
    bld.program(source=['src/pvm2sdp/main.cxx',
                        'src/pvm2sdp/parse_command_line.cxx',
                        'src/pvm2sdp/read_input_files/read_input_files.cxx',
                        'src/pvm2sdp/read_input_files/Input_Parser/on_start_element.cxx',
                        'src/pvm2sdp/read_input_files/Input_Parser/on_end_element.cxx',
                        'src/pvm2sdp/read_input_files/Input_Parser/on_characters.cxx',
                        'src/pvm2sdp/Dual_Constraint_Group/Dual_Constraint_Group/Dual_Constraint_Group.cxx',
                        'src/pvm2sdp/Dual_Constraint_Group/Dual_Constraint_Group/sample_bilinear_basis.cxx',
                        'src/pvm2sdp/write_objectives.cxx',
                        'src/pvm2sdp/write_bilinear_bases.cxx',
                        'src/pvm2sdp/write_blocks.cxx',
                        'src/pvm2sdp/write_primal_objective_c.cxx',
                        'src/pvm2sdp/write_free_var_matrix.cxx'],
                target='pvm2sdp',
                cxxflags=default_flags,
                use=use_packages
                )

    bld.program(source=['src/sdp2blocks/main.cxx',
                        'src/sdp2blocks/read_input/read_input.cxx',
                        'src/sdp2blocks/read_input/Input_Parser/on_start_element.cxx',
                        'src/sdp2blocks/read_input/Input_Parser/on_end_element.cxx',
                        'src/sdp2blocks/read_input/Input_Parser/on_characters.cxx',
                        'src/sdp2blocks/write_output_files/write_output_files.cxx',
                        'src/sdp2blocks/write_output_files/sample_points.cxx',
                        'src/sdp2blocks/write_output_files/bilinear_basis/bilinear_basis.cxx',
                        'src/sdp2blocks/write_output_files/bilinear_basis/bilinear_form/bilinear_form.cxx',
                        'src/sdp2blocks/write_output_files/bilinear_basis/bilinear_form/integral.cxx',
                        'src/sdp2blocks/write_output_files/bilinear_basis/bilinear_form/rest.cxx',
                        'src/sdp2blocks/write_output_files/bilinear_basis/bilinear_form/dExp.cxx',
                        'src/sdp2blocks/write_output_files/bilinear_basis/bilinear_form/derivative.cxx',
                        'src/sdp2blocks/write_output_files/bilinear_basis/bilinear_form/operator_plus_set_Derivative_Term.cxx'],
                target='sdp2blocks',
                cxxflags=default_flags,
                use=use_packages
                )

    bld.program(source=['src/block_grid_mapping/main.cxx',
                        'src/compute_block_grid_mapping.cxx'],
                target='block_grid_mapping',
                cxxflags=default_flags,
                use=use_packages
                )
                
                        
    
