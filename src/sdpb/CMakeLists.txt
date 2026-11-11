add_executable(sdpb
  main.cxx
  SDPB_Parameters.cxx)
target_link_libraries(sdpb PRIVATE sdp_solve sdpa_solve MPSolve::MPSolve)
