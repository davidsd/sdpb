# GPU backend: the device library and its checks, with the sources and flags of src/gpu/Makefile.
# Built with -DSDPB_ENABLE_GPU=ON, see CMakeLists.txt in the root directory.

# Options of every CUDA compilation here (NVFLAGS in the Makefile).
# -UNDEBUG: keep assert() enabled, as in the Makefile and waf builds.
# -fPIC, as in waf builds: the C++ compiler that links the objects may make position-independent executables.
add_library(sdpb_gpu_options INTERFACE)
target_compile_features(sdpb_gpu_options INTERFACE cuda_std_17)
target_compile_options(sdpb_gpu_options INTERFACE "$<$<COMPILE_LANGUAGE:CUDA>:-Xcompiler=-Wall;-Xcompiler=-fopenmp;-Xcompiler=-fPIC;-UNDEBUG>")
target_link_libraries(sdpb_gpu_options INTERFACE GMP::gmp)

# Libraries of the device code (LDLIBS in the Makefile)
add_library(sdpb_gpu_libs INTERFACE)
target_link_libraries(sdpb_gpu_libs INTERFACE
  CUDA::cudart_static
  CUDA::cublasLt
  GMP::gmp
  CBLAS::CBLAS
  $<LINK_ONLY:OpenMP::OpenMP_CXX>)

add_library(sdpb_gpu STATIC
  rns/Rns_Basis.cxx
  rns/Gram_Accumulator.cu
  rns/Mpf_Gram.cu
  rns/Rns_Gemm.cu
  linalg/Block_Factor.cu
  solver/Mpf_Ops.cu
  solver/Gpu_Solver.cu)
# compiled by nvcc too, as in the Makefile and waf builds
set_source_files_properties(rns/Rns_Basis.cxx PROPERTIES LANGUAGE CUDA)
target_link_libraries(sdpb_gpu PRIVATE sdpb_gpu_options PUBLIC sdpb_gpu_libs)

# Checks (src/gpu/README.md, "Tests"), outside the default build, in this build directory, e.g.
#   cmake --build build --target gpu_checks mpfr_checks

# Self-contained or on dumps of real blocks; they need a GPU
add_executable(gram_check EXCLUDE_FROM_ALL test/gram_check.cu)
target_link_libraries(gram_check PRIVATE sdpb_gpu sdpb_gpu_options)
add_executable(mpf_gram_check EXCLUDE_FROM_ALL test/mpf_gram_check.cu)
target_link_libraries(mpf_gram_check PRIVATE sdpb_gpu sdpb_gpu_options)
add_executable(gemm_check EXCLUDE_FROM_ALL test/gemm_check.cu)
target_link_libraries(gemm_check PRIVATE sdpb_gpu sdpb_gpu_options)
add_executable(chol_tile_check EXCLUDE_FROM_ALL test/chol_tile_check.cu)
target_link_libraries(chol_tile_check PRIVATE sdpb_gpu_options sdpb_gpu_libs)
add_executable(mpf_gemm_check EXCLUDE_FROM_ALL test/mpf_gemm_check.cu)
target_link_libraries(mpf_gemm_check PRIVATE sdpb_gpu sdpb_gpu_options)
add_executable(block_factor_check EXCLUDE_FROM_ALL test/block_factor_check.cu)
target_link_libraries(block_factor_check PRIVATE sdpb_gpu sdpb_gpu_options)
add_executable(real_block_check EXCLUDE_FROM_ALL test/real_block_check.cu)
target_link_libraries(real_block_check PRIVATE sdpb_gpu sdpb_gpu_options)
add_executable(q_check EXCLUDE_FROM_ALL test/q_check.cu)
target_link_libraries(q_check PRIVATE sdpb_gpu sdpb_gpu_options)
add_executable(residual_check EXCLUDE_FROM_ALL test/residual_check.cu)
target_link_libraries(residual_check PRIVATE sdpb_gpu sdpb_gpu_options)
add_executable(mpf_ops_check EXCLUDE_FROM_ALL test/mpf_ops_check.cu)
target_link_libraries(mpf_ops_check PRIVATE sdpb_gpu sdpb_gpu_options)
add_executable(factor_profile EXCLUDE_FROM_ALL test/factor_profile.cu)
target_link_libraries(factor_profile PRIVATE sdpb_gpu sdpb_gpu_options)
add_executable(gemv_compare EXCLUDE_FROM_ALL test/gemv_compare.cu)
target_link_libraries(gemv_compare PRIVATE sdpb_gpu sdpb_gpu_options)
add_custom_target(gpu_checks DEPENDS gram_check mpf_gram_check gemm_check chol_tile_check mpf_gemm_check block_factor_check real_block_check q_check residual_check mpf_ops_check factor_profile gemv_compare)

# MPFR references for the exported solves, on the CPU
add_executable(ref_solve EXCLUDE_FROM_ALL test/ref_solve.c)
target_compile_options(ref_solve PRIVATE -UNDEBUG)
target_link_libraries(ref_solve PRIVATE OpenMP::OpenMP_C MPFR::MPFR GMP::gmp m)
add_executable(compare_references EXCLUDE_FROM_ALL test/compare_references.c)
target_compile_options(compare_references PRIVATE -UNDEBUG)
target_link_libraries(compare_references PRIVATE MPFR::MPFR GMP::gmp m)
add_custom_target(mpfr_checks DEPENDS ref_solve compare_references)

set_target_properties(gram_check mpf_gram_check gemm_check chol_tile_check mpf_gemm_check block_factor_check real_block_check q_check residual_check mpf_ops_check factor_profile gemv_compare PROPERTIES
  RUNTIME_OUTPUT_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}"
  # linked by the C++ compiler, as sdpb
  LINKER_LANGUAGE CXX)
set_target_properties(ref_solve compare_references PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}")
