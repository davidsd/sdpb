add_library(sdp_solve STATIC
  Solver_Parameters/Memory_Limit.cxx
  Solver_Parameters/Solver_Parameters.cxx
  Block_Info/Block_Info.cxx
  Block_Matrix/Block_Diagonal_Matrix.cxx
  SDP/SDP.cxx
  SDP/read_normalization.cxx
  SDP/read_objectives.cxx
  SDP/assign_bilinear_bases_dist.cxx
  SDP/read_block_data/read_block_data.cxx
  SDP/read_block_data/SDP_Block_Data.cxx
  SDP/set_bases_blocks.cxx
  SDP_Solver/save_checkpoint.cxx
  SDP_Solver/load_checkpoint/load_checkpoint.cxx
  SDP_Solver/SDP_Solver.cxx
  SDP_Solver/run/run.cxx
  SDP_Solver/run/constraint_matrix_weighted_sum.cxx
  SDP_Solver/run/compute_dual_residues_and_error.cxx
  SDP_Solver/run/compute_primal_residues_and_error_P_Ax_X.cxx
  SDP_Solver/run/compute_primal_residues_and_error_p_b_Bx.cxx
  SDP_Solver/run/compute_objectives.cxx
  SDP_Solver/run/compute_bilinear_pairings/compute_bilinear_pairings.cxx
  SDP_Solver/run/compute_bilinear_pairings/compute_A_X_inv.cxx
  SDP_Solver/run/compute_bilinear_pairings/compute_A_Y.cxx
  SDP_Solver/run/compute_feasible_and_termination.cxx
  SDP_Solver/run/print_header.cxx
  SDP_Solver/run/print_iteration.cxx
  SDP_Solver/run/step/step.cxx
  SDP_Solver/run/step/initialize_schur_complement_solver/initialize_schur_complement_solver.cxx
  SDP_Solver/run/step/initialize_schur_complement_solver/compute_schur_complement.cxx
  SDP_Solver/run/step/initialize_schur_complement_solver/compute_Q.cxx
  SDP_Solver/run/bigint_syrk/BigInt_Shared_Memory_Syrk_Context/BigInt_Shared_Memory_Syrk_Context.cxx
  SDP_Solver/run/bigint_syrk/BigInt_Shared_Memory_Syrk_Context/bigint_syrk_blas.cxx
  SDP_Solver/run/bigint_syrk/BigInt_Shared_Memory_Syrk_Context/clear_residues.cxx
  SDP_Solver/run/bigint_syrk/BigInt_Shared_Memory_Syrk_Context/compute_block_residues.cxx
  SDP_Solver/run/bigint_syrk/BigInt_Shared_Memory_Syrk_Context/get_blas_job_schedule.cxx
  SDP_Solver/run/bigint_syrk/BigInt_Shared_Memory_Syrk_Context/restore_and_reduce.cxx
  SDP_Solver/run/bigint_syrk/blas_jobs/Blas_Job.cxx
  SDP_Solver/run/bigint_syrk/blas_jobs/create_blas_job_schedule.cxx
  SDP_Solver/run/bigint_syrk/Bigint_Syrk_Config.cxx
  SDP_Solver/run/step/compute_search_direction/compute_search_direction.cxx
  SDP_Solver/run/step/compute_search_direction/compute_schur_RHS.cxx
  SDP_Solver/run/step/compute_search_direction/solve_schur_complement_equation.cxx
  SDP_Solver_Terminate_Reason.cxx
  Write_Solution.cxx)
target_link_libraries(sdp_solve PUBLIC sdpb_util pmp2sdp_lib)

if(SDPB_ENABLE_GPU)
  target_sources(sdp_solve PRIVATE
    SDP_Solver/run/gpu/Gpu_Compute.cxx)
  target_link_libraries(sdp_solve PUBLIC sdpb_gpu)
endif()
