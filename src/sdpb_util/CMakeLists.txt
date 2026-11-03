add_library(sdpb_util STATIC
  Archive_Reader.cxx
  Boost_Float.cxx
  copy_matrix.cxx
  Environment.cxx
  memory_estimates.cxx
  Memory_Tracker.cxx
  Mesh.cxx
  Proc_Meminfo.cxx
  bigint_shared_memory/blas_jobs/Blas_Job_Cost.cxx
  bigint_shared_memory/fmpz/Fmpz_BigInt.cxx
  bigint_shared_memory/fmpz/Fmpz_Matrix.cxx
  bigint_shared_memory/fmpz/Fmpz_Comb.cxx
  bigint_shared_memory/fmpz/residues.cxx
  block_mapping/Abstract_Block_Info.cxx
  block_mapping/Block_Mapping.cxx
  Timers/Scoped_Timer.cxx
  Timers/Timer.cxx
  Timers/Timers.cxx)
target_link_libraries(sdpb_util PUBLIC sdpb_options ${SDPB_DEPENDENCIES})
