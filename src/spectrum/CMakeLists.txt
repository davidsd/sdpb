add_executable(spectrum
  main.cxx
  handle_arguments.cxx
  read_c_minus_By.cxx
  read_pmp_info.cxx
  read_x.cxx
  set_default_parameters.cxx
  write_profiling.cxx
  compute_spectrum/compute_spectrum.cxx
  compute_spectrum/find_zeros.cxx
  compute_spectrum/mpsolve.cxx
  write_spectrum/write_spectrum.cxx
  write_spectrum/write_file.cxx)
target_link_libraries(spectrum PRIVATE sdpb_util MPSolve::MPSolve)
