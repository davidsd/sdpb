add_executable(approx_objective
  main.cxx
  Approx_Parameters.cxx
  Axpy.cxx
  setup_solver.cxx
  write_solver_state.cxx
  Approx_Objective/Approx_Objective.cxx
  Approx_Objective/compute_dx_dy.cxx
  linear_approximate_objectives.cxx
  quadratic_approximate_objectives.cxx)
target_link_libraries(approx_objective PRIVATE pmp_read sdp_solve)
