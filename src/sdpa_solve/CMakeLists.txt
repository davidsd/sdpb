add_library(sdpa_solve STATIC
  SDP.cxx
  Block_Info.cxx
  parse_sdpa/parse_sdpa.cxx
  SDP_Solver/save_checkpoint.cxx
  SDP_Solver/run/print_header.cxx
  SDP_Solver/run/compute_dual_residues_and_error.cxx
  SDP_Solver/run/step/compute_S.cxx
  SDP_Solver/run/step/Solver_Run_Config.cxx
  SDP_Solver/run/step/compute_search_direction.cxx
  SDP_Solver/run/step/Initialize_P_Config.cxx
  SDP_Solver/run/step/Initialize_P_Context.cxx
  SDP_Solver/run/step/step.cxx
  SDP_Solver/run/step/bigint_trmm/blas_jobs/Blas_Job.cxx
  SDP_Solver/run/step/bigint_trmm/blas_jobs/create_blas_job_schedule.cxx
  SDP_Solver/run/run.cxx
  SDP_Solver/run/compute_primal_residues_and_error.cxx
  SDP_Solver/run/print_iteration.cxx
  SDP_Solver/run/compute_objectives.cxx
  SDP_Solver/load_checkpoint/load_checkpoint.cxx
  SDP_Solver/SDP_Solver.cxx)
target_link_libraries(sdpa_solve PUBLIC sdp_solve)
