add_library(catch2 STATIC "${PROJECT_SOURCE_DIR}/external/catch2/catch_amalgamated.cpp")
target_include_directories(catch2 PUBLIC "${PROJECT_SOURCE_DIR}/external")
target_compile_definitions(catch2 PUBLIC CATCH_AMALGAMATED_CUSTOM_MAIN)

add_executable(integration_tests
  src/integration_tests/main.cxx
  src/integration_tests/util/Float.cxx
  src/integration_tests/util/diff_outer_limits.cxx
  src/integration_tests/util/diff_sdp.cxx
  src/integration_tests/util/diff_sdpb_out.cxx
  src/integration_tests/util/diff_spectrum.cxx
  src/integration_tests/util/Parse_SDP.cxx
  src/integration_tests/util/process.cxx
  src/integration_tests/util/Test_Case_Runner.cxx
  src/integration_tests/cases/end-to-end.test.cxx
  src/integration_tests/cases/outer_limits.test.cxx
  src/integration_tests/cases/pmp2sdp.test.cxx
  src/integration_tests/cases/SDPLIB.test.cxx
  src/integration_tests/cases/resident.test.cxx
  src/integration_tests/cases/sdpb.test.cxx)
target_include_directories(integration_tests PRIVATE src)
target_link_libraries(integration_tests PRIVATE sdpb_util catch2)
if(TARGET Boost::process)
  target_link_libraries(integration_tests PRIVATE Boost::process)
endif()

add_executable(unit_tests
  src/unit_tests/main.cxx
  src/unit_tests/cases/LPT_scheduling.test.cxx
  src/unit_tests/cases/Matrix_Normalizer.test.cxx
  src/unit_tests/cases/Memory_Tracker.test.cxx
  src/unit_tests/cases/block_data_serialization.test.cxx
  src/unit_tests/cases/block_mapping.test.cxx
  src/unit_tests/cases/Boost_Float.test.cxx
  src/unit_tests/cases/blas_residues.test.cxx
  src/unit_tests/cases/boost_serialization.test.cxx
  src/unit_tests/cases/create_blas_job_schedule.test.cxx
  src/unit_tests/cases/calculate_matrix_square.test.cxx
  src/unit_tests/cases/copy_matrix.test.cxx
  src/unit_tests/cases/json.test.cxx
  src/unit_tests/cases/multimodular_blas.test.cxx
  src/unit_tests/cases/pmp_sampling.test.cxx
  src/unit_tests/cases/shared_window.test.cxx
  src/unit_tests/cases/Sparse_Block_Matrix.test.cxx)
target_include_directories(unit_tests PRIVATE src)
target_link_libraries(unit_tests PRIVATE pmp_read pmp2sdp_lib sdp_solve catch2)

# Tests for ctest

# MPI launcher, must accept "-n <processes>"
set(SDPB_TEST_MPIRUN "" CACHE STRING "MPI launcher for tests, e.g. 'srun --mpi=pmi2'. Empty: use mpiexec")
set(mpirun_command "${SDPB_TEST_MPIRUN}")
if(NOT mpirun_command)
  if(MPIEXEC_EXECUTABLE)
    set(mpirun_command "${MPIEXEC_EXECUTABLE}")
    execute_process(COMMAND "${MPIEXEC_EXECUTABLE}" --version OUTPUT_VARIABLE version ERROR_QUIET)
    # Open MPI refuses to start more processes than CPU cores
    if(version MATCHES "Open MPI|OpenRTE")
      string(APPEND mpirun_command " --oversubscribe")
    endif()
  else()
    # Not found by CMake, look for it in PATH when running tests
    set(mpirun_command mpiexec)
  endif()
endif()
message(STATUS "MPI launcher for tests: ${mpirun_command}")

separate_arguments(mpirun UNIX_COMMAND "${mpirun_command}")
# calculate_matrix_square test requires a multiple of 6 processes
add_test(NAME unit_tests COMMAND ${mpirun} -n 6 $<TARGET_FILE:unit_tests>)
# Reads test/data in the source directory, writes test/out in the build directory
add_test(NAME integration_tests COMMAND integration_tests "--mpirun=${mpirun_command}"
  WORKING_DIRECTORY "${PROJECT_SOURCE_DIR}")
