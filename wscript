import os

def options(opt):
    opt.load(['compiler_cxx','gnu_dirs','boost','gmpxx','cxx14','elemental',
              'libxmlxx'])

def configure(conf):
    if not 'CXX' in os.environ or os.environ['CXX']=='g++' or os.environ['CXX']=='icpc':
        conf.environ['CXX']='mpicxx'

    conf.load(['compiler_cxx','gnu_dirs','boost','gmpxx','cxx14','elemental',
               'libxmlxx'])
    conf.check_boost(lib='serialization system filesystem timer program_options chrono')

def build(bld):
    default_flags=['-Wall', '-Wextra', '-O3', '-Wno-deprecated']
    # default_flags=['-Wall', '-Wextra', '-g', '-Wno-deprecated']
    use_packages=['BOOST','gmpxx','cxx14','elemental','libxmlxx']
    
    # Main executable
    bld.program(source=['src/main.cxx',
                        'src/SDP_Solver_Parameters/ostream.cxx',
                        'src/solve/solve.cxx',
                        'src/solve/SDP_Solver/save_checkpoint.cxx',
                        'src/solve/SDP_Solver/load_checkpoint.cxx',
                        'src/solve/SDP_Solver/SDP_Solver.cxx',
                        'src/solve/SDP/SDP/SDP.cxx',
                        'src/solve/SDP/SDP/Input_Parser/on_start_element.cxx',
                        'src/solve/SDP/SDP/Input_Parser/on_end_element.cxx',
                        'src/solve/SDP/SDP/Input_Parser/on_characters.cxx',
                        'src/solve/SDP/SDP/parse_vector.cxx',
                        'src/solve/SDP/SDP/parse_BigFloat.cxx',
                        'src/solve/SDP/SDP/parse_polynomial_vector_matrix.cxx',
                        'src/solve/SDP/SDP/bootstrap/bootstrap.cxx',
                        'src/solve/SDP/SDP/bootstrap/dual_constraint_group_from_pol_vec_mat/dual_constraint_group_from_pol_vec_mat.cxx',
                        'src/solve/SDP/SDP/bootstrap/dual_constraint_group_from_pol_vec_mat/sample_bilinear_basis.cxx',
                        'src/solve/SDP/SDP/bootstrap/fill_from_dual_constraint_groups/fill_from_dual_constraint_groups.cxx',
                        'src/solve/SDP/SDP/bootstrap/fill_from_dual_constraint_groups/set_block_sizes.cxx',
                        'src/solve/SDP/SDP/bootstrap/fill_from_dual_constraint_groups/set_bilinear.cxx',
                        'src/solve/SDP/SDP/bootstrap/fill_from_dual_constraint_groups/assign_blocks.cxx',
                        'src/solve/SDP/SDP/bootstrap/fill_from_dual_constraint_groups/compute_block_grid_mapping.cxx',
                        'src/solve/SDP_Solver/run/constraint_matrix_weighted_sum.cxx',
                        'src/solve/SDP_Solver/run/compute_dual_residues.cxx',
                        'src/solve/SDP_Solver/run/dot.cxx',
                        'src/solve/SDP_Solver/run/block_tensor_inv_transpose_congruence_with_cholesky.cxx',
                        'src/solve/SDP_Solver/run/corrector_centering_parameter.cxx',
                        'src/solve/SDP_Solver/run/predictor_centering_parameter.cxx',
                        'src/solve/SDP_Solver/run/run.cxx',
                        'src/solve/SDP_Solver/run/step_length.cxx',
                        'src/solve/SDP_Solver/run/compute_primal_residues.cxx',
                        'src/solve/SDP_Solver/run/block_tensor_transpose_congruence.cxx',
                        'src/solve/SDP_Solver/initialize_schur_complement_solver/initialize_schur_complement_solver.cxx',
                        'src/solve/SDP_Solver/initialize_schur_complement_solver/compute_schur_complement.cxx',
                        'src/solve/SDP_Solver/print_iteration.cxx',
                        'src/solve/SDP_Solver/save_solution.cxx',
                        'src/solve/SDP_Solver/run/compute_search_direction/compute_search_direction.cxx',
                        'src/solve/SDP_Solver/run/compute_search_direction/compute_schur_RHS.cxx',
                        'src/solve/SDP_Solver/run/compute_search_direction/solve_schur_complement_equation.cxx',
                        'src/solve/SDP_Solver/print_header.cxx',
                        'src/solve/SDP_Solver_Terminate_Reason/ostream.cxx',
                        'src/solve/lower_triangular_solve.cxx',
                        'src/solve/lower_triangular_transpose_solve.cxx',
                        'src/solve/Block_Diagonal_Matrix/block_diagonal_matrix_multiply.cxx',
                        'src/solve/Block_Diagonal_Matrix/block_diagonal_matrix_scale_multiply_add.cxx',
                        'src/solve/Block_Diagonal_Matrix/block_matrix_solve_with_cholesky.cxx',
                        'src/solve/Block_Diagonal_Matrix/cholesky_decomposition.cxx',
                        'src/solve/Block_Diagonal_Matrix/frobenius_product_of_sums.cxx',
                        'src/solve/Block_Diagonal_Matrix/frobenius_product_symmetric.cxx',
                        'src/solve/Block_Diagonal_Matrix/lower_triangular_inverse_congruence.cxx',
                        'src/solve/Block_Diagonal_Matrix/min_eigenvalue.cxx',
                        'src/solve/Block_Diagonal_Matrix/ostream.cxx'],
                target='sdpb',
                cxxflags=default_flags,
                rpath=[bld.env.LIBDIR],
                use=use_packages
                )
