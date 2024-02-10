# TODO choose version:
cmake_minimum_required(VERSION 3.16)
project(sdpb LANGUAGES CXX)

# This is to ensure that CMake also checks lib64 for the dependencies 
set(CMAKE_PREFIX_PATH "/usr/lib64" ${CMAKE_PREFIX_PATH})

# Currently we are working with a single CMake file but it is better to refactor parts of this file and put them inside cmake/modules
list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")
list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake/modules")


# cxx 17
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
#set(CMAKE_CXX_EXTENSIONS OFF)

# We all love some colour don't we
include(CustomMessages)

# compiler_cxx
# TODO
if(DEFINED (ENV{CXX}))
    if(ENV{CXX} == "g++" || ENV{CXX} == "icpc")
        SET(CMAKE_CXX_COMPILER mpicxx)
    endif()
else()
    SET(CMAKE_CXX_COMPILER mpicxx)
endif()

message(STATUS "Found g++ compiler (C++ compiler): ${CMAKE_CXX_COMPILER}")

# TODO gnu_dirs

##################################################
########### Find required dependencies ###########
##################################################

include(FindPkgConfig)
find_package(PkgConfig REQUIRED)

# TODO boost
# TODO Check that the custom approach works. Have not tested.

find_package(Boost REQUIRED COMPONENTS system filesystem date_time program_options iostreams serialization)
if(NOT Boost_FOUND)
        pkg_check_modules(Boost REQUIRED COMPONENTS system filesystem date_time program_options iostreams serialization)
endif()
message(STATUS "Boost version: ${Boost_VERSION}")

# Omitted the part where boost.py tests for boost by running a small snippet of code.

# TODO gmpxx
# NOTE This does not work right now because find_package is not able to find FindGMPXX.cmake
# find_package(GMPXX REQUIRED)
# message(STATUS "GMP version: ${GMPXX_VERSION}")


# TODO mpfr
find_package(MPFR QUIET)

if (MPFR_FOUND) 
	message(STATUS "MPFR version: ${MPFR_VERSION}")
else ()
	find_path(MPFR_INCLUDE_DIR mpfr.h
		PATHS /usr/include /usr/local/include
	)

	find_library(MPFR_LIBRARY mpfr
    		PATHS /usr/lib64 /usr/local/lib
	)

	if (MPFR_INCLUDE_DIR AND MPFR_LIBRARY)
    		set(MPFR_FOUND TRUE)
	endif ()

	if (MPFR_FOUND)
    		if (NOT MPFR_FIND_QUIETLY)
        		message(STATUS "Found MPFR: ${MPFR_LIBRARY}")
    		endif ()
	else ()
    		if (MPFR_FIND_REQUIRED)
        		message(FATAL_ERROR "Could not find MPFR")
    		endif ()
	endif ()

	mark_as_advanced(MPFR_INCLUDE_DIR MPFR_LIBRARY)

	message(STATUS "MPFR version: ${MPFR_VERSION}")
endif ()

# TODO elemental


# TODO libxml2
#find_package(libxml2 REQUIRED)
# if(NOT libxml2_FOUND)
#         pkg_check_modules(LIBXML2 REQUIRED libxml-2)
# endif()
# message(STATUS "libxml2 version: ${libxml2_VERSION}")

# TODO rapidjson
pkg_check_modules(RapidJSON REQUIRED RapidJSON)
message(STATUS "RapidJSON version: ${RapidJSON_VERSION}")

# TODO libarchive
pkg_check_modules(libarchive REQUIRED libarchive)
message(STATUS "libarchive version: ${libarchive_VERSION}")


#def configure(conf):
#    if not 'CXX' in os.environ or os.environ['CXX']=='g++' or os.environ['CXX']=='icpc':
#        conf.environ['CXX']='mpicxx'

find_package(MPI REQUIRED)
message(STATUS "Run: ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} ${MPIEXEC_MAX_NUMPROCS} ${MPIEXEC_PREFLAGS} EXECUTABLE ${MPIEXEC_POSTFLAGS} ARGS")

# TODO Refactor into a separate file for Git (see Elemental for example)
#    conf.env.git_version=subprocess.check_output('git describe --dirty', universal_newlines=True, shell=True).rstrip()
find_package(Git QUIET REQUIRED)
if(NOT EXISTS "${PROJECT_SOURCE_DIR}/.git")
    message(FATAL_ERROR "${PROJECT_SOURCE_DIR}/.git not found")
endif()
execute_process(COMMAND ${GIT_EXECUTABLE} describe --always --dirty
        WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
        RESULT_VARIABLE GIT_DESCRIBE_RESULT
        OUTPUT_VARIABLE SDPB_VERSION_STRING
        OUTPUT_STRIP_TRAILING_WHITESPACE)
if(NOT GIT_DESCRIBE_RESULT EQUAL "0")
    message(FATAL_ERROR "git describe --always --dirty failed: ${GIT_DESCRIBE_RESULT}.")
endif()

#def build(bld):
#    default_flags=['-Wall', '-Wextra', '-O3', '-DOMPI_SKIP_MPICXX', '-D SDPB_VERSION_STRING="' + bld.env.git_version + '"']
#    # default_flags=['-Wall', '-Wextra', '-O3', '-g', '-DOMPI_SKIP_MPICXX', '-D SDPB_VERSION_STRING="' + bld.env.git_version + '"']
#    # default_flags=['-Wall', '-Wextra', '-g', '-DOMPI_SKIP_MPICXX', '-D SDPB_VERSION_STRING="' + bld.env.git_version + '"']

#TODO -O3 etc is set by choosing cmake debug/release
# TODO -Werror=return-type, but it requires to mark El::RuntimeError as [[noreturn]]
add_compile_options(-Wall -Wextra -DOMPI_SKIP_MPICXX)
add_compile_definitions(SDPB_VERSION_STRING="${SDPB_VERSION_STRING}")
add_compile_definitions(RAPIDJSON_HAS_STDSTRING=1)

# TODO setup flags for build configurations (debug/release)

if(CMAKE_BUILD_TYPE STREQUAL "Release")
  # This option is okay as-is
  set(CMAKE_BUILD_TYPE Release)
elseif(CMAKE_BUILD_TYPE STREQUAL "Debug")
  # This option is okay as-is
  set(CMAKE_BUILD_TYPE Debug)
else()
  message(WARNING "Build mode not specified, defaulting to Release build.")
  set(CMAKE_BUILD_TYPE Release)
endif()

if(CMAKE_BUILD_TYPE STREQUAL "Release")
  set(SDPB_RELEASE TRUE)
else()
  set(SDPB_RELEASE FALSE)
endif()
string(TOUPPER ${CMAKE_BUILD_TYPE} UPPER_BUILD_TYPE)


#    use_packages=['cxx17','boost','gmpxx','mpfr','elemental','libxml2', 'rapidjson', 'libarchive']
#TODO

set(CMAKE_LIBRARY_PATH CMAKE_LIBRARY_PATH /usr/local/lib)
message("CMAKE_LIBRARY_PATH = ${CMAKE_LIBRARY_PATH}")

include(FeatureSummary)

set(libxml2_INCLUDE_DIR /usr/include/libxml2)
include(FindLibXml2)
#find_package(libxml2)
#find_library(LIBXML2 libxml2 REQUIRED)

find_library(MPFR_LIBRARIES mpfr)

#find_library()
# TODO use target_link_libraries
#link_libraries(cxx17 boost gmpxx mpfr elemental libxml2 rapidjson libarchive)

find_library(LIBARCHIVE_LIBRARIES archive)

# TODO This does not work right now and gives undefined variable error, because of which cmake cannot configure the build
find_library(elemental elemental)
link_libraries(elemental)
foreach(LIB El pmrrr ElSuiteSparse)
    find_library(Elemental_${LIB}_LIBRARIES ${LIB}
            PATHS $ENV{HOME}/install
            PATH_SUFFIXES lib lib64)
    list(APPEND ELEMENTAL_LIBRARIES "${Elemental_${LIB}_LIBRARIES}")
endforeach()

# see elemental/cmake/modules/FindGMP.cmake
#find_path(GMP_INCLUDES NAMES gmp.h PATHS $ENV{GMPDIR} $ENV{GMP_HOME} $ENV{GMPHOME} $ENV{GMP_INCLUDE} ${INCLUDE_INSTALL_DIR})
find_library(GMP_LIBRARIES gmp PATHS $ENV{GMPHOME} $ENV{GMPDIR} $ENV{GMP_HOME} $ENV{GMP_LIB} ${LIB_INSTALL_DIR})

# see elemental/cmake/modules/FindGMPXX.cmake
#find_path(GMPXX_INCLUDES NAMES gmpxx.h PATHS $ENV{GMPXXDIR} $ENV{GMPXX_HOME} $ENV{GMPXX_INCLUDE} $ENV{GMPDIR} $ENV{GMP_HOME} $ENV{GMP_INCLUDE} ${INCLUDE_INSTALL_DIR})
find_library(GMPXX_LIBRARIES gmpxx PATHS $ENV{GMPHOME} $ENV{GMPXXDIR} $ENV{GMPXX_HOME} $ENV{GMPXX_LIB} $ENV{GMPDIR} $ENV{GMP_HOME} $ENV{GMP_LIB} ${LIB_INSTALL_DIR})

# FIXME This is a hardcoded path
# NOTE FLINT is not required in the current release of SDPB. Commenting it out...
#find_library(FLINT_LIBRARIES flint PATHS /home/vasdommes/install/lib)

# From FLINT - FindCBLAS
find_path(CBLAS_INCLUDE_DIRS NAMES cblas.h
        HINTS CBLAS_ROOT ENV CBLAS_ROOT
        PATHS ${INCLUDE_INSTALL_DIR} ${CMAKE_INSTALL_PREFIX}/include
        PATH_SUFFIXES openblas cblas blis
)

find_library(CBLAS_LIBRARIES NAMES accelerate openblas cblas blas blis
        HINTS CBLAS_ROOT ENV CBLAS_ROOT
        PATHS ${LIB_INSTALL_DIR} ${CMAKE_INSTALL_PREFIX}/lib
        PATH_SUFFIXES openblas cblas blis
)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(CBLAS
        "Could NOT find a BLAS compatible library or 'cblas.h', install BLAS or set CBLAS_ROOT."
        CBLAS_INCLUDE_DIRS CBLAS_LIBRARIES)

mark_as_advanced(CBLAS_LIBRARIES CBLAS_INCLUDE_DIRS)


include_directories(src)
include_directories(external)
include_directories(SYSTEM elemental/install/include)
include_directories(SYSTEM ${libxml2_INCLUDE_DIR})
include_directories(SYSTEM /home/vasdommes/install/include)

add_library(sdp_solve)
target_sources(sdp_solve
        #TODO set correct visibility
        PUBLIC
        src/sdp_solve/Solver_Parameters/Solver_Parameters.cxx
        src/sdp_solve/Solver_Parameters/ostream.cxx
        src/sdp_solve/Solver_Parameters/to_property_tree.cxx
        src/sdp_solve/Archive_Reader/Archive_Reader.cxx
        src/sdp_solve/Archive_Reader/underflow.cxx
        src/sdp_solve/Block_Info/Block_Info.cxx
        src/sdp_solve/Block_Info/read_block_info.cxx
        src/sdp_solve/Block_Info/read_block_costs.cxx
        src/sdp_solve/Block_Info/allocate_blocks/allocate_blocks.cxx
        src/sdp_solve/Block_Info/allocate_blocks/compute_block_grid_mapping.cxx
        src/sdp_solve/SDP/SDP/SDP.cxx
        src/sdp_solve/SDP/SDP/read_objectives.cxx
        src/sdp_solve/SDP/SDP/assign_bilinear_bases_dist.cxx
        src/sdp_solve/SDP/SDP/read_blocks/read_blocks.cxx
        src/sdp_solve/SDP/SDP/read_blocks/read_block_stream/read_block_stream.cxx
        src/sdp_solve/SDP/SDP/read_blocks/read_block_stream/Block_Parser/EndArray.cxx
        src/sdp_solve/SDP/SDP/read_blocks/read_block_stream/Block_Parser/Key.cxx
        src/sdp_solve/SDP/SDP/read_blocks/read_block_stream/Block_Parser/StartArray.cxx
        src/sdp_solve/SDP/SDP/read_blocks/read_block_stream/Block_Parser/String.cxx
        src/sdp_solve/SDP/SDP/set_bases_blocks.cxx
        src/sdp_solve/SDP_Solver/save_checkpoint.cxx
        src/sdp_solve/SDP_Solver/load_checkpoint/load_checkpoint.cxx
        src/sdp_solve/SDP_Solver/load_checkpoint/load_binary_checkpoint.cxx
        src/sdp_solve/SDP_Solver/load_checkpoint/load_text_checkpoint.cxx
        src/sdp_solve/SDP_Solver/SDP_Solver.cxx
        src/sdp_solve/SDP_Solver/run/run.cxx
        src/sdp_solve/SDP_Solver/run/cholesky_decomposition.cxx
        src/sdp_solve/SDP_Solver/run/constraint_matrix_weighted_sum.cxx
        src/sdp_solve/SDP_Solver/run/compute_dual_residues_and_error.cxx
        src/sdp_solve/SDP_Solver/run/compute_primal_residues_and_error_P_Ax_X.cxx
        src/sdp_solve/SDP_Solver/run/compute_primal_residues_and_error_p_b_Bx.cxx
        src/sdp_solve/SDP_Solver/run/compute_objectives/compute_objectives.cxx
        src/sdp_solve/SDP_Solver/run/compute_objectives/dot.cxx
        src/sdp_solve/SDP_Solver/run/compute_bilinear_pairings/compute_bilinear_pairings.cxx
        src/sdp_solve/SDP_Solver/run/compute_bilinear_pairings/compute_A_X_inv.cxx
        src/sdp_solve/SDP_Solver/run/compute_bilinear_pairings/compute_A_Y.cxx
        src/sdp_solve/SDP_Solver/run/compute_feasible_and_termination.cxx
        src/sdp_solve/SDP_Solver/run/print_header.cxx
        src/sdp_solve/SDP_Solver/run/print_iteration.cxx
        src/sdp_solve/SDP_Solver/run/step/step.cxx
        src/sdp_solve/SDP_Solver/run/step/initialize_schur_complement_solver/initialize_schur_complement_solver.cxx
        src/sdp_solve/SDP_Solver/run/step/initialize_schur_complement_solver/compute_schur_complement.cxx
        #src/sdp_solve/SDP_Solver/run/step/initialize_schur_complement_solver/compute_Q.cxx
        #src/sdp_solve/SDP_Solver/run/bigint_syrk/BigInt_Shared_Memory_Syrk_Context.cxx
        #src/sdp_solve/SDP_Solver/run/bigint_syrk/compute_block_residues.cxx
        #src/sdp_solve/SDP_Solver/run/bigint_syrk/Matrix_Normalizer.cxx
        #src/sdp_solve/SDP_Solver/run/bigint_syrk/Fmpz_Matrix.cxx
        #src/sdp_solve/SDP_Solver/run/bigint_syrk/Fmpz_Comb.cxx
        #src/sdp_solve/SDP_Solver/run/bigint_syrk/blas_jobs/Blas_Job.cxx
        #src/sdp_solve/SDP_Solver/run/bigint_syrk/blas_jobs/Blas_Job_Schedule.cxx
        #src/sdp_solve/SDP_Solver/run/bigint_syrk/blas_jobs/create_blas_job_schedule.cxx
        src/sdp_solve/SDP_Solver/run/step/compute_search_direction/compute_search_direction.cxx
        src/sdp_solve/SDP_Solver/run/step/compute_search_direction/cholesky_solve.cxx
        src/sdp_solve/SDP_Solver/run/step/compute_search_direction/compute_schur_RHS.cxx
        src/sdp_solve/SDP_Solver/run/step/compute_search_direction/scale_multiply_add.cxx
        src/sdp_solve/SDP_Solver/run/step/compute_search_direction/solve_schur_complement_equation.cxx
        src/sdp_solve/SDP_Solver/run/step/predictor_centering_parameter.cxx
        src/sdp_solve/SDP_Solver/run/step/corrector_centering_parameter/corrector_centering_parameter.cxx
        src/sdp_solve/SDP_Solver/run/step/corrector_centering_parameter/frobenius_product_of_sums.cxx
        src/sdp_solve/SDP_Solver/run/step/frobenius_product_symmetric.cxx
        src/sdp_solve/SDP_Solver/run/step/step_length/step_length.cxx
        src/sdp_solve/SDP_Solver/run/step/step_length/min_eigenvalue.cxx
        src/sdp_solve/SDP_Solver/run/step/step_length/lower_triangular_inverse_congruence.cxx
        src/sdp_solve/SDP_Solver_Terminate_Reason/ostream.cxx
        src/sdp_solve/lower_triangular_transpose_solve.cxx
        src/sdp_solve/Block_Diagonal_Matrix/ostream.cxx
        src/sdp_solve/Write_Solution.cxx)
# TODO add_subdirectory
target_link_libraries(sdp_solve
        Boost::program_options
        Boost::serialization
        ${MPFR_LIBRARIES}
        ${LIBARCHIVE_LIBRARIES}
        ${GMP_LIBRARIES} ${GMPXX_LIBRARIES}
        ${ELEMENTAL_LIBRARIES}
        #${FLINT_LIBRARIES}
        ${CBLAS_LIBRARIES})


# SDPB executable
add_executable(sdpb)
target_sources(sdpb
        PUBLIC
        src/sdpb/main.cxx
        src/sdpb/solve.cxx
        src/sdpb/write_timing.cxx
        src/sdpb/SDPB_Parameters/SDPB_Parameters.cxx
        src/sdpb/SDPB_Parameters/to_property_tree.cxx
        src/sdpb/SDPB_Parameters/ostream.cxx
        src/sdpb/save_solution.cxx)
target_link_libraries(sdpb
        sdp_solve
        Boost::program_options
        Boost::serialization
        ${MPFR_LIBRARIES}
        ${GMP_LIBRARIES} ${GMPXX_LIBRARIES}
        ${ELEMENTAL_LIBRARIES})

add_library(sdp_convert)
target_sources(sdp_convert
        PUBLIC
        src/sdp_convert/Dual_Constraint_Group/Dual_Constraint_Group.cxx
        src/sdp_convert/Dual_Constraint_Group/sample_bilinear_basis.cxx
        src/sdp_convert/write_block_data.cxx
        src/sdp_convert/write_block_info_json.cxx
        src/sdp_convert/write_objectives_json.cxx
        src/sdp_convert/write_sdpb_input_files.cxx
        src/sdp_convert/write_control_json.cxx
        src/sdp_convert/Archive_Writer/Archive_Writer.cxx
        src/sdp_convert/Archive_Writer/write_entry.cxx
        src/sdp_convert/Archive_Entry.cxx)
target_link_libraries(sdp_convert
        Boost::iostreams
        Boost::serialization
        ${LIBARCHIVE_LIBRARIES})

add_library(sdp_read)
target_sources(sdp_read
        PUBLIC
        src/sdp_read/read_input/read_input.cxx
        src/sdp_read/read_input/read_json/read_json.cxx
        src/sdp_read/read_input/read_json/Positive_Matrix_With_Prefactor_State/json_key.cxx
        src/sdp_read/read_input/read_json/Positive_Matrix_With_Prefactor_State/json_string.cxx
        src/sdp_read/read_input/read_json/Positive_Matrix_With_Prefactor_State/json_start_array.cxx
        src/sdp_read/read_input/read_json/Positive_Matrix_With_Prefactor_State/json_end_array.cxx
        src/sdp_read/read_input/read_json/Positive_Matrix_With_Prefactor_State/json_start_object.cxx
        src/sdp_read/read_input/read_json/Positive_Matrix_With_Prefactor_State/json_end_object.cxx
        src/sdp_read/read_input/read_json/Damped_Rational_State/json_key.cxx
        src/sdp_read/read_input/read_json/Damped_Rational_State/json_string.cxx
        src/sdp_read/read_input/read_json/Damped_Rational_State/json_start_array.cxx
        src/sdp_read/read_input/read_json/Damped_Rational_State/json_end_array.cxx
        src/sdp_read/read_input/read_json/Damped_Rational_State/json_start_object.cxx
        src/sdp_read/read_input/read_json/Damped_Rational_State/json_end_object.cxx
        src/sdp_read/read_input/read_json/JSON_Parser/Key.cxx
        src/sdp_read/read_input/read_json/JSON_Parser/String.cxx
        src/sdp_read/read_input/read_json/JSON_Parser/StartArray.cxx
        src/sdp_read/read_input/read_json/JSON_Parser/EndArray.cxx
        src/sdp_read/read_input/read_json/JSON_Parser/StartObject.cxx
        src/sdp_read/read_input/read_json/JSON_Parser/EndObject.cxx
        src/sdp_read/read_input/read_mathematica/read_mathematica.cxx
        src/sdp_read/read_input/read_mathematica/parse_SDP/parse_SDP.cxx
        src/sdp_read/read_input/read_mathematica/parse_SDP/parse_matrices.cxx
        src/sdp_read/read_input/read_mathematica/parse_SDP/parse_number.cxx
        src/sdp_read/read_input/read_mathematica/parse_SDP/parse_polynomial.cxx
        src/sdp_read/read_input/read_mathematica/parse_SDP/parse_matrix/parse_matrix.cxx
        src/sdp_read/read_input/read_mathematica/parse_SDP/parse_matrix/parse_damped_rational.cxx
        src/sdp_read/read_pvm_input/read_pvm_input.cxx
        src/sdp_read/read_pvm_input/read_xml_input/read_xml_input.cxx
        src/sdp_read/read_pvm_input/read_xml_input/Input_Parser/on_start_element.cxx
        src/sdp_read/read_pvm_input/read_xml_input/Input_Parser/on_end_element.cxx
        src/sdp_read/read_pvm_input/read_xml_input/Input_Parser/on_characters.cxx
        src/sdp_read/read_nsv_file_list.cxx
        src/sdp_read/sample_points.cxx
        src/sdp_read/sample_scalings.cxx)
target_link_libraries(sdp_read sdp_convert)

add_executable(pvm2sdp)
target_sources(pvm2sdp
        PUBLIC
        src/pvm2sdp/main.cxx
        src/pvm2sdp/parse_command_line.cxx
        src/pvm2sdp/read_input_files/read_input_files.cxx
        src/pvm2sdp/read_input_files/read_xml_input/read_xml_input.cxx
        src/pvm2sdp/read_input_files/read_xml_input/Input_Parser/on_start_element.cxx
        src/pvm2sdp/read_input_files/read_xml_input/Input_Parser/on_end_element.cxx
        src/pvm2sdp/read_input_files/read_xml_input/Input_Parser/on_characters.cxx)
target_link_libraries(pvm2sdp sdp_read)

add_executable(sdp2input)
target_sources(sdp2input
        PUBLIC
        src/sdp2input/main.cxx
        src/sdp2input/write_output/write_output.cxx
        src/sdp2input/write_output/bilinear_basis/bilinear_basis.cxx
        src/sdp2input/write_output/bilinear_basis/precompute/precompute.cxx
        src/sdp2input/write_output/bilinear_basis/precompute/integral.cxx
        src/sdp2input/write_output/bilinear_basis/bilinear_form/bilinear_form.cxx
        src/sdp2input/write_output/bilinear_basis/bilinear_form/rest.cxx
        src/sdp2input/write_output/bilinear_basis/bilinear_form/dExp.cxx
        src/sdp2input/write_output/bilinear_basis/bilinear_form/derivative.cxx
        src/sdp2input/write_output/bilinear_basis/bilinear_form/operator_plus_set_Derivative_Term.cxx)
target_link_libraries(sdp2input sdp_read)

add_library(mesh)
target_sources(mesh
        PUBLIC
        src/Mesh/Mesh.cxx
        src/Mesh/ostream.cxx)

add_executable(outer_limits)
target_sources(outer_limits
        PUBLIC
        src/outer_limits/main.cxx
        src/outer_limits/power_prefactor.cxx
        src/outer_limits/poles_prefactor.cxx
        src/outer_limits/Function/eval/eval.cxx
        src/outer_limits/compute_optimal/compute_optimal.cxx
        src/outer_limits/compute_optimal/compute_y_transform.cxx
        src/outer_limits/compute_optimal/setup_constraints.cxx
        src/outer_limits/compute_optimal/copy_matrix.cxx
        src/outer_limits/compute_optimal/find_new_points/find_new_points.cxx
        src/outer_limits/compute_optimal/find_new_points/eval_summed.cxx
        src/outer_limits/compute_optimal/find_new_points/get_new_points.cxx
        src/outer_limits/compute_optimal/load_checkpoint/load_checkpoint.cxx
        src/outer_limits/compute_optimal/load_checkpoint/Checkpoint_Parser/EndArray.cxx
        src/outer_limits/compute_optimal/load_checkpoint/Checkpoint_Parser/EndObject.cxx
        src/outer_limits/compute_optimal/load_checkpoint/Checkpoint_Parser/Key.cxx
        src/outer_limits/compute_optimal/load_checkpoint/Checkpoint_Parser/StartArray.cxx
        src/outer_limits/compute_optimal/load_checkpoint/Checkpoint_Parser/StartObject.cxx
        src/outer_limits/compute_optimal/load_checkpoint/Checkpoint_Parser/String.cxx
        src/outer_limits/compute_optimal/save_checkpoint.cxx
        src/outer_limits/read_points/read_points.cxx
        src/outer_limits/read_points/read_points_json/read_points_json.cxx
        src/outer_limits/read_points/read_points_json/Points_Parser/EndArray.cxx
        src/outer_limits/read_points/read_points_json/Points_Parser/EndObject.cxx
        src/outer_limits/read_points/read_points_json/Points_Parser/Key.cxx
        src/outer_limits/read_points/read_points_json/Points_Parser/StartArray.cxx
        src/outer_limits/read_points/read_points_json/Points_Parser/StartObject.cxx
        src/outer_limits/read_points/read_points_json/Points_Parser/String.cxx
        src/outer_limits/read_function_blocks/read_function_blocks.cxx
        src/outer_limits/read_function_blocks/read_json/Function_State/json_end_array.cxx
        src/outer_limits/read_function_blocks/read_json/Function_State/json_end_object.cxx
        src/outer_limits/read_function_blocks/read_json/Function_State/json_key.cxx
        src/outer_limits/read_function_blocks/read_json/Function_State/json_start_array.cxx
        src/outer_limits/read_function_blocks/read_json/Function_State/json_start_object.cxx
        src/outer_limits/read_function_blocks/read_json/Function_State/json_string.cxx
        src/outer_limits/read_function_blocks/read_json/Function_Blocks_Parser/EndArray.cxx
        src/outer_limits/read_function_blocks/read_json/Function_Blocks_Parser/EndObject.cxx
        src/outer_limits/read_function_blocks/read_json/Function_Blocks_Parser/Key.cxx
        src/outer_limits/read_function_blocks/read_json/Function_Blocks_Parser/StartArray.cxx
        src/outer_limits/read_function_blocks/read_json/Function_Blocks_Parser/StartObject.cxx
        src/outer_limits/read_function_blocks/read_json/Function_Blocks_Parser/String.cxx
        src/outer_limits/read_function_blocks/read_json/read_json.cxx
        src/outer_limits/Outer_Parameters/Outer_Parameters.cxx
        src/outer_limits/Outer_Parameters/to_property_tree.cxx
        src/outer_limits/Outer_Parameters/ostream.cxx)
target_link_libraries(outer_limits sdp_read sdp_solve mesh)

add_executable(approx_objective)
target_sources(approx_objective
        PUBLIC
        src/approx_objective/main.cxx
        src/approx_objective/Approx_Parameters/Approx_Parameters.cxx
        src/approx_objective/Approx_Parameters/ostream.cxx
        src/approx_objective/Axpy.cxx
        src/approx_objective/setup_solver.cxx
        src/approx_objective/write_solver_state.cxx
        src/approx_objective/Approx_Objective/Approx_Objective/Approx_Objective.cxx
        src/approx_objective/Approx_Objective/Approx_Objective/compute_dx_dy.cxx
        src/approx_objective/linear_approximate_objectives.cxx
        src/approx_objective/quadratic_approximate_objectives.cxx)
target_link_libraries(approx_objective sdp_read sdp_solve)

add_executable(pvm2functions)
target_sources(pvm2functions
        PUBLIC
        src/pvm2functions/main.cxx
        src/pvm2functions/parse_command_line.cxx
        src/pvm2functions/write_functions.cxx)
target_link_libraries(pvm2functions sdp_read)

add_executable(sdp2functions)
target_sources(sdp2functions
        PUBLIC
        src/sdp2functions/main.cxx
        src/sdp2functions/write_functions.cxx)
target_link_libraries(sdp2functions sdp_read)

add_executable(spectrum)
target_sources(spectrum
        PUBLIC
        src/spectrum/main.cxx
        src/spectrum/handle_arguments.cxx
        src/spectrum/read_x.cxx
        src/spectrum/compute_spectrum_pmp.cxx
        src/spectrum/compute_spectrum_pvm.cxx
        src/spectrum/compute_lambda.cxx
        src/spectrum/eval_summed.cxx
        src/spectrum/get_zeros.cxx
        src/spectrum/write_spectrum/write_spectrum.cxx
        src/spectrum/write_spectrum/write_file.cxx)
target_link_libraries(spectrum sdp_read sdp_solve sdp_convert mesh)

add_executable(unit_tests)
target_compile_definitions(unit_tests PUBLIC CATCH_AMALGAMATED_CUSTOM_MAIN)
target_include_directories(unit_tests PUBLIC src test/src)
target_sources(unit_tests
        PUBLIC
        external/catch2/catch_amalgamated.cpp
        test/src/unit_tests/main.cxx
        test/src/unit_tests/cases/block_data_serialization.test.cxx
        test/src/unit_tests/cases/boost_serialization.test.cxx
        #test/src/unit_tests/cases/calculate_matrix_square.test.cxx
        #test/src/unit_tests/cases/create_blas_job_schedule.test.cxx
        #test/src/unit_tests/cases/LPT_scheduling.test.cxx
        #test/src/unit_tests/cases/Matrix_Normalizer.test.cxx
        #test/src/unit_tests/cases/shared_window.test.cxx
        )
target_link_libraries(unit_tests
        sdp_convert
        sdp_solve
        ${MPFR_LIBRARIES}
        ${GMP_LIBRARIES} ${GMPXX_LIBRARIES}
        ${ELEMENTAL_LIBRARIES})

add_executable(integration_tests)
target_compile_definitions(integration_tests PUBLIC CATCH_AMALGAMATED_CUSTOM_MAIN)
target_include_directories(integration_tests PUBLIC src test/src SYSTEM ${GMPXX_INCLUDES})
target_sources(integration_tests
        PUBLIC
        external/catch2/catch_amalgamated.cpp
        test/src/integration_tests/main.cxx
        test/src/integration_tests/util/Float.cxx
        test/src/integration_tests/util/diff_outer_limits.cxx
        test/src/integration_tests/util/diff_sdp_zip.cxx
        test/src/integration_tests/util/diff_sdpb_out.cxx
        test/src/integration_tests/util/diff_spectrum.cxx
        test/src/integration_tests/util/Test_Case_Runner.cxx
        test/src/integration_tests/cases/end-to-end.test.cxx
        test/src/integration_tests/cases/outer_limits.test.cxx
        test/src/integration_tests/cases/pvm2sdp.test.cxx
        test/src/integration_tests/cases/sdp2input.test.cxx
        test/src/integration_tests/cases/sdpb.test.cxx
        test/src/integration_tests/cases/spectrum.test.cxx)
target_link_libraries(integration_tests
        Boost::serialization
        Boost::filesystem
        ${MPFR_LIBRARIES}
        ${GMP_LIBRARIES} ${GMPXX_LIBRARIES}
        ${ELEMENTAL_LIBRARIES}
)
