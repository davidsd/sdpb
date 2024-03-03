# TODO choose version, currently working with 3.16 as a temporary measure:
cmake_minimum_required(VERSION 3.16)

project(sdpb LANGUAGES CXX)
# TODO setting version of SDPB

# TODO - Currently we are working with a single CMake file but it is better to refactor parts of this file and put them inside cmake/modules
list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")
list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake/modules")

# Custom colour scheme for logging
include(CustomMessages)

# This is to ensure that CMake also checks lib64 for the dependencies on Linux systems
if(UNIX)
        set(CMAKE_PREFIX_PATH "/usr/lib64" ${CMAKE_PREFIX_PATH})
endif()

# cxx 17
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
#set(CMAKE_CXX_EXTENSIONS OFF)

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

# GMP and GMPXX
# -------------
find_package(GMP QUIET)

if (GMP_FOUND) 
	message(STATUS "GMP version: ${GMP_VERSION}")
else ()
	find_path(GMP_INCLUDE_DIR 
                        NAMES gmp.h 
                        PATHS $ENV{GMPDIR} $ENV{GMP_HOME} $ENV{GMP_INCLUDE} ${INCLUDE_INSTALL_DIR}
                        )

	find_library(GMP_LIBRARIES gmp 
                        PATHS $ENV{GMPDIR} $ENV{GMP_HOME} $ENV{GMP_LIB} ${LIB_INSTALL_DIR}
                        )

	if (GMP_INCLUDE_DIR AND GMP_LIBRARIES)
    		set(GMP_FOUND TRUE)
	endif ()

	if (GMP_FOUND)
    		if (NOT GMP_FIND_QUIETLY)
        		message(STATUS "Found GMP: ${GMP_LIBRARIES}")
    		endif ()
	else ()
    		if (GMP_FIND_REQUIRED)
        		message(FATAL_ERROR "Could not find GMP")
    		endif ()
	endif ()

	mark_as_advanced(GMP_INCLUDE_DIR GMP_LIBRARIES)

	message(STATUS "GMP version: ${GMP_VERSION}")
endif ()
if(GMP_FOUND AND GMP_INCLUDE_DIR)
    include_directories(${GMP_INCLUDE_DIR})
endif()

find_package(GMPXX QUIET)

if (GMPXX_FOUND) 
	message(STATUS "GMPXX version: ${GMPXX_VERSION}")
else ()
	find_path(GMPXX_INCLUDE_DIR 
                        NAMES gmpxx.h 
                        PATHS $ENV{GMPXXDIR} $ENV{GMPXX_HOME} $ENV{GMPXX_INCLUDE} $ENV{GMPDIR} $ENV{GMP_HOME} $ENV{GMP_INCLUDE} ${INCLUDE_INSTALL_DIR}
                        )

	find_library(GMPXX_LIBRARIES gmpxx
                        PATHS $ENV{GMPXXDIR} $ENV{GMPXX_HOME} $ENV{GMPXX_LIB} $ENV{GMPDIR} $ENV{GMP_HOME} $ENV{GMP_LIB} ${LIB_INSTALL_DIR}
                        )

	if (GMPXX_INCLUDE_DIR AND GMPXX_LIBRARIES)
    		set(GMPXX_FOUND TRUE)
	endif ()

	if (GMPXX_FOUND)
    		if (NOT GMPXX_FIND_QUIETLY)
        		message(STATUS "Found GMPXX: ${GMPXX_LIBRARIES}")
    		endif ()
	else ()
    		if (GMPXX_FIND_REQUIRED)
        		message(FATAL_ERROR "Could not find GMPXX")
    		endif ()
	endif ()

	mark_as_advanced(GMPXX_INCLUDE_DIR GMPXX_LIBRARIES)

	message(STATUS "GMPXX version: ${GMPXX_VERSION}")
endif ()
if(GMPXX_FOUND AND GMPXX_INCLUDE_DIR)
    include_directories(${GMPXX_INCLUDES})
endif()

# MPFR
# ----
find_package(MPFR QUIET)

if (MPFR_FOUND) 
	message(STATUS "MPFR version: ${MPFR_VERSION}")
else ()
	find_path(MPFR_INCLUDE_DIR mpfr.h
		PATHS /usr/include /usr/local/include
	)

	find_library(MPFR_LIBRARIES mpfr
    		PATHS /usr/lib64 /usr/local/lib
	)

	if (MPFR_INCLUDE_DIR AND MPFR_LIBRARIES)
    		set(MPFR_FOUND TRUE)
	endif ()

	if (MPFR_FOUND)
    		if (NOT MPFR_FIND_QUIETLY)
        		message(STATUS "Found MPFR: ${MPFR_LIBRARIES}")
    		endif ()
	else ()
    		if (MPFR_FIND_REQUIRED)
        		message(FATAL_ERROR "Could not find MPFR")
    		endif ()
	endif ()

	mark_as_advanced(MPFR_INCLUDE_DIR MPFR_LIBRARIES)

	message(STATUS "MPFR version: ${MPFR_VERSION}")
endif ()
message(STATUS ${MPFR_INCLUDE_DIR})

if(MPFR_FOUND AND MPFR_INCLUDE_DIR)
    include_directories(${MPFR_INCLUDE_DIR})
endif()


# Boost
# -----
# TODO Change pkg_check_modules to a solution where you use find_path and find_library 

find_package(Boost 
                REQUIRED COMPONENTS 
                system 
                filesystem 
                date_time 
                program_options 
                iostreams 
                serialization 
                stacktrace_basic
                stacktrace_backtrace
                stacktrace_addr2line
                stacktrace_noop)
if(NOT Boost_FOUND)
        pkg_check_modules(Boost REQUIRED COMPONENTS system filesystem date_time program_options iostreams serialization)
endif()
message(STATUS "Boost version: ${Boost_VERSION}")
if(Boost_FOUND AND Boost_INCLUDE_DIR)
    include_directories(${Boost_INCLUDE_DIR})
endif()

# Omitted the part where boost.py tests for boost by running a small snippet of code.
# TODO - Maybe add the tests back again?

# libxml2
# -------
# TODO Change pkg_check_modules to a solution where you use find_path and find_library
find_package(LibXml2 REQUIRED)
if(NOT LIBXML2_FOUND)
        pkg_check_modules(LIBXML2 REQUIRED libxml-2)
endif()
message(STATUS "libxml2 version: ${LIBXML2_VERSION} at ${LIBXML2_DIRECTORY}")
if(LIBXML2_FOUND AND LIBXML2_INCLUDE_DIR)
    include_directories(${LIBXML2_INCLUDE_DIR})
endif()


# RapidJSON
# ---------
# TODO Change pkg_check_modules to a solution where you use find_path and find_library
find_package(RapidJSON REQUIRED)
# pkg_check_modules(RapidJSON REQUIRED RapidJSON)
message(STATUS "RapidJSON version: ${RapidJSON_VERSION}")
if(RapidJSON_FOUND AND RapidJSON_INCLUDE_DIR)
    include_directories(${RapidJSON_INCLUDE_DIR})
endif()


# libarchive
# ----------
# TODO Change pkg_check_modules to a solution where you use find_path and find_library
find_package(LibArchive REQUIRED)
if (NOT libarchive_FOUND)
        pkg_check_modules(libarchive REQUIRED libarchive)
endif()      
message(STATUS "libarchive version: ${libarchive_VERSION}")
if(LibArchive_FOUND AND LibArchive_INCLUDE_DIR)
    include_directories(${LibArchive_INCLUDE_DIR})
endif()


#def configure(conf):
#    if not 'CXX' in os.environ or os.environ['CXX']=='g++' or os.environ['CXX']=='icpc':
#        conf.environ['CXX']='mpicxx'

# MPI
# ---
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

# set(CMAKE_LIBRARY_PATH CMAKE_LIBRARY_PATH /usr/local/lib)
# message("CMAKE_LIBRARY_PATH = ${CMAKE_LIBRARY_PATH}")

include(FeatureSummary)

#find_library()
# TODO use target_link_libraries
#link_libraries(cxx17 boost gmpxx mpfr elemental libxml2 rapidjson libarchive)

# Elemental
# ---------
find_package(Elemental REQUIRED)
# find_path(Elemental_INCLUDES El.hpp
#         PATHS /usr/include /usr/local/include ${CMAKE_PREFIX_PATH}
# )

# find_library(Elemental_LIBRARIES El
#         PATHS /usr/lib64 /usr/local/lib ${CMAKE_PREFIX_PATH}
# )

# if (Elemental_INCLUDES AND Elemental_LIBRARIES)
#         set(Elemental_FOUND TRUE)
# endif ()

# if (Elemental_FOUND)
#         if (NOT Elemental_FIND_QUIETLY)
#                 message(STATUS "Found El: ${Elemental_LIBRARIES}")
#         endif ()
# else ()
#         if (Elemental_FIND_REQUIRED)
#                 message(FATAL_ERROR "Could not find Elemental")
#         endif ()
# endif ()
# mark_as_advanced(Elemental_INCLUDES Elemental_LIBRARIES)

if(Elemental_FOUND AND Elemental_INCLUDE_DIRS)
    include_directories(${Elemental_INCLUDE_DIRS})
endif()

message(STATUS "El version: ${Elemental_VERSION}")
message(STATUS "Found Elemental: ${Elemental_LIBRARIES}")
# find_library(El El PATHS ${Elemental_DIR} REQUIRED)
find_library(PMRRR_LIBRARIES pmrrr PATHS ${Elemental_DIR}/external/pmrrr REQUIRED)
find_library(ELSUITESPARSE_LIBRARIES ElSuiteSparse PATHS ${Elemental_DIR}/external/suite_sparse REQUIRED)
list(APPEND Elemental_LIBRARIES ${PMRRR_LIBRARIES})
list(APPEND Elemental_LIBRARIES ${ELSUITESPARSE_LIBRARIES})
message(STATUS "Elemental Libraries - ${Elemental_LIBRARIES}")

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

set(USE_PACKAGES ${GMP_LIBRARIES} ${GMPXX_LIBRARIES} ${MPFR_LIBRARIES} ${Elemental_LIBRARIES} ${LIBXML2_LIBRARIES} ${LibArchive_LIBRARIES})
message(STATUS ${USE_PACKAGES})

add_subdirectory(src/sdpb_util)

add_subdirectory(src/pmp)

add_subdirectory(src/pmp2sdp)

add_subdirectory(src/pmp_read)

add_subdirectory(src/sdp2input)

add_subdirectory(src/pvm2sdp)

# Currently difficult to put the add_executable call to src/pmp2sdp/CMakeLists.txt because between pmp2sdp_lib and pmp2sdp, pmp_read needs to be built
# TODO - Refactor the code below once a resolution to the above conflict has been figured out

add_executable(pmp2sdp)
target_sources(pmp2sdp
                PRIVATE
                src/pmp2sdp/main.cxx
                src/pmp2sdp/Pmp2sdp_Parameters/Pmp2sdp_Parameters.cxx)
target_link_libraries(pmp2sdp 
                        PRIVATE
                        ${USE_PACKAGES} 
                        Boost::iostreams
                        Boost::serialization
                        Boost::stacktrace_basic
                        Boost::stacktrace_addr2line
                        Boost::stacktrace_noop
                        Boost::program_options
                        dl
                        sdpb_util
                        pmp 
                        pmp_read 
                        pmp2sdp_lib)

add_subdirectory(src/pmp2functions)

add_subdirectory(src/sdp_solve)

add_subdirectory(src/sdpb)

add_subdirectory(src/outer_limits)

add_subdirectory(src/approx_objective)

add_subdirectory(src/spectrum)

# add_executable(integration_tests)
# target_sources(integration_tests
#                 PRIVATE
#                 external/catch2/catch_amalgamated.cpp
#                 test/src/integration_tests/main.cxx
#                 test/src/integration_tests/util/Float.cxx
#                 test/src/integration_tests/util/diff_outer_limits.cxx
#                 test/src/integration_tests/util/diff_sdp.cxx
#                 test/src/integration_tests/util/diff_sdpb_out.cxx
#                 test/src/integration_tests/util/diff_spectrum.cxx
#                 test/src/integration_tests/util/Test_Case_Runner.cxx
#                 test/src/integration_tests/cases/end-to-end.test.cxx
#                 test/src/integration_tests/cases/outer_limits.test.cxx
#                 test/src/integration_tests/cases/pmp2sdp.test.cxx
#                 test/src/integration_tests/cases/sdpb.test.cxx)
# target_link_libraries(integration_tests
#                         PRIVATE
#                         ${USE_PACKAGES}
#                         Boost::serialization)
# target_include_directories(integration_tests PRIVATE test/src)
#     bld.program(source=['external/catch2/catch_amalgamated.cpp',
#                         'test/src/integration_tests/main.cxx',
#                         'test/src/integration_tests/util/Float.cxx',
#                         'test/src/integration_tests/util/diff_outer_limits.cxx',
#                         'test/src/integration_tests/util/diff_sdp.cxx',
#                         'test/src/integration_tests/util/diff_sdpb_out.cxx',
#                         'test/src/integration_tests/util/diff_spectrum.cxx',
#                         'test/src/integration_tests/util/Test_Case_Runner.cxx',
#                         'test/src/integration_tests/cases/end-to-end.test.cxx',
#                         'test/src/integration_tests/cases/outer_limits.test.cxx',
#                         'test/src/integration_tests/cases/pmp2sdp.test.cxx',
#                         'test/src/integration_tests/cases/sdpb.test.cxx'],
#                 target='integration_tests',
#                 install_path=None,
#                 cxxflags=default_flags,
#                 defines=default_defines + ['CATCH_AMALGAMATED_CUSTOM_MAIN'],
#                 use=use_packages,
#                 includes=default_includes + ['test/src']
#                 )
#     bld.program(source=['external/catch2/catch_amalgamated.cpp',
#                         'test/src/unit_tests/main.cxx',
#                         'test/src/unit_tests/cases/block_data_serialization.test.cxx',
#                         'test/src/unit_tests/cases/block_mapping.test.cxx',
#                         'test/src/unit_tests/cases/Boost_Float.test.cxx',
#                         'test/src/unit_tests/cases/boost_serialization.test.cxx',
#                         'test/src/unit_tests/cases/copy_matrix.test.cxx',
#                         'test/src/unit_tests/cases/json.test.cxx',
#                         'test/src/unit_tests/cases/shared_window.test.cxx'],
#                 target='unit_tests',
#                 cxxflags=default_flags,
#                 defines=default_defines + ['CATCH_AMALGAMATED_CUSTOM_MAIN'],
#                 use=use_packages + ['pmp_read', 'pmp2sdp_lib', 'sdp_solve'],
#                 includes=default_includes + ['test/src']
#                 )