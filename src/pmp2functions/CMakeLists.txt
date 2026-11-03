add_executable(pmp2functions
  main.cxx
  Pmp2functions_Parameters.cxx
  write_functions.cxx)
target_link_libraries(pmp2functions PRIVATE pmp_read)
