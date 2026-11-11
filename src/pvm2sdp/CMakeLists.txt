add_executable(pvm2sdp
  main.cxx
  parse_command_line.cxx)
target_link_libraries(pvm2sdp PRIVATE pmp_read)
