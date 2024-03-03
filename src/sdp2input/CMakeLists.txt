add_executable(sdp2input)
target_sources(sdp2input PRIVATE main.cxx)
target_link_libraries(sdp2input 
                        PUBLIC
                        ${USE_PACKAGES} 
                        Boost::program_options
                        Boost::iostreams
                        Boost::serialization
                        Boost::stacktrace_basic
                        Boost::stacktrace_addr2line
                        Boost::stacktrace_noop
                        dl
                        pmp 
                        pmp_read 
                        pmp2sdp_lib
                        sdpb_util)