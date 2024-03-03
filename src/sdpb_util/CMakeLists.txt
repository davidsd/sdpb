add_library(sdpb_util SHARED)
target_sources(sdpb_util
                PRIVATE
                copy_matrix.cxx
                Environment.cxx
                Mesh.cxx
                Proc_Meminfo.cxx
                Timers/Scoped_Timer.cxx
                Timers/Timer.cxx
                Timers/Timers.cxx)
target_link_libraries(sdpb_util PUBLIC ${GMPXX_LIBRARIES} Boost::iostreams Boost::serialization Boost::filesystem ${Elemental_LIBRARIES})
