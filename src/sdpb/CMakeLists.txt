add_executable(sdpb
  main.cxx
  SDPB_Parameters.cxx
  resident/Control_Channel.cxx
  resident/Resident.cxx)
target_link_libraries(sdpb PRIVATE sdp_solve sdpa_solve pmp_read MPSolve::MPSolve)
