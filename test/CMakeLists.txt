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
  src/unit_tests/cases/block_data_serialization.test.cxx
  src/unit_tests/cases/block_mapping.test.cxx
  src/unit_tests/cases/Boost_Float.test.cxx
  src/unit_tests/cases/boost_serialization.test.cxx
  src/unit_tests/cases/create_blas_job_schedule.test.cxx
  src/unit_tests/cases/calculate_matrix_square.test.cxx
  src/unit_tests/cases/copy_matrix.test.cxx
  src/unit_tests/cases/json.test.cxx
  src/unit_tests/cases/pmp_sampling.test.cxx
  src/unit_tests/cases/shared_window.test.cxx)
target_include_directories(unit_tests PRIVATE src)
target_link_libraries(unit_tests PRIVATE pmp_read pmp2sdp_lib sdp_solve catch2)
