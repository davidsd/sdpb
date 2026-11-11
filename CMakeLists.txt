cmake_minimum_required(VERSION 3.21)
project(SDPB LANGUAGES CXX)

if(PROJECT_SOURCE_DIR STREQUAL PROJECT_BINARY_DIR)
  message(FATAL_ERROR "In-source builds are not supported. Use a build directory, e.g. 'cmake -B build'.")
endif()
if(EXISTS "${PROJECT_BINARY_DIR}/c4che")
  message(FATAL_ERROR "${PROJECT_BINARY_DIR} contains a waf build. Remove it or use another build directory.")
endif()

list(APPEND CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake")

if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type: Release, Debug, RelWithDebInfo or MinSizeRel" FORCE)
endif()
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
# All executables go to the build directory, e.g. build/sdpb
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${PROJECT_BINARY_DIR}")
# Installed executables find dependencies outside system directories
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH ON)

# Dependencies

# SDPB does not use the deprecated MPI C++ bindings
set(MPI_CXX_SKIP_MPICXX ON)
find_package(MPI REQUIRED COMPONENTS CXX)
find_package(Boost 1.70 REQUIRED CONFIG
  COMPONENTS date_time filesystem iostreams program_options serialization
  OPTIONAL_COMPONENTS stacktrace_backtrace stacktrace_addr2line stacktrace_basic)
# Boost.Process, used by integration tests, is a compiled library since Boost 1.89
if(Boost_VERSION VERSION_GREATER_EQUAL 1.89)
  find_package(Boost REQUIRED CONFIG COMPONENTS process)
endif()
find_package(GMP REQUIRED)
find_package(MPFR REQUIRED)
find_package(Elemental REQUIRED)
find_package(LibXml2 REQUIRED)
find_package(RapidJSON REQUIRED)
find_package(LibArchive REQUIRED)
find_package(FLINT 2.8 REQUIRED)
find_package(CBLAS REQUIRED)
find_package(MPSolve 3.2.2 REQUIRED)

# Boost.Stacktrace prints source locations only with a backend library.
# Use the first backend that links, otherwise header-only.
include(CheckCXXSourceCompiles)
# Check results are cached, repeat the checks if Boost changes
if(NOT Boost_DIR STREQUAL "${_SDPB_STACKTRACE_CHECKED_BOOST_DIR}")
  foreach(backend IN ITEMS backtrace addr2line basic)
    unset(SDPB_BOOST_STACKTRACE_${backend}_LINKS CACHE)
  endforeach()
  set(_SDPB_STACKTRACE_CHECKED_BOOST_DIR "${Boost_DIR}" CACHE INTERNAL "")
endif()
set(stacktrace_backend ${CMAKE_DL_LIBS})
set(stacktrace_define "")
foreach(backend IN ITEMS backtrace addr2line basic)
  if(NOT Boost_stacktrace_${backend}_FOUND)
    continue()
  endif()
  set(CMAKE_REQUIRED_LIBRARIES Boost::stacktrace_${backend} ${CMAKE_DL_LIBS})
  set(CMAKE_REQUIRED_DEFINITIONS -DBOOST_STACKTRACE_LINK)
  check_cxx_source_compiles([[
#include <boost/stacktrace.hpp>
int main() { boost::stacktrace::stacktrace(); }
]] SDPB_BOOST_STACKTRACE_${backend}_LINKS)
  unset(CMAKE_REQUIRED_LIBRARIES)
  unset(CMAKE_REQUIRED_DEFINITIONS)
  if(SDPB_BOOST_STACKTRACE_${backend}_LINKS)
    set(stacktrace_backend Boost::stacktrace_${backend} ${CMAKE_DL_LIBS})
    set(stacktrace_define BOOST_STACKTRACE_LINK)
    break()
  endif()
endforeach()
message(STATUS "Boost.Stacktrace backend: ${stacktrace_backend}")

# External dependencies of all SDPB targets except MPSolve.
# Include directories are searched in this order.
set(SDPB_DEPENDENCIES
  MPI::MPI_CXX
  GMP::gmpxx
  MPFR::MPFR
  Boost::headers Boost::date_time Boost::filesystem Boost::iostreams Boost::program_options Boost::serialization
  ${stacktrace_backend}
  Elemental::El
  LibXml2::LibXml2
  RapidJSON::RapidJSON
  LibArchive::LibArchive
  FLINT::FLINT
  CBLAS::CBLAS)

include(CheckPackageVersions)
sdpb_check_package_versions(${SDPB_DEPENDENCIES} MPSolve::MPSolve)

# Generates sdpb_version.hxx at every build
set(SDPB_VERSION_STRING "" CACHE STRING "SDPB version, e.g. 3.1.0. Empty: run 'git describe' at every build")
find_package(Git QUIET)
set(generated_dir "${PROJECT_BINARY_DIR}/generated")
add_custom_target(sdpb_version
  COMMAND "${CMAKE_COMMAND}" "-DVERSION=${SDPB_VERSION_STRING}"
    -DSOURCE_DIR=${PROJECT_SOURCE_DIR} -DGIT_EXECUTABLE=${GIT_EXECUTABLE}
    -DOUTPUT=${generated_dir}/sdpb_version.hxx -P ${PROJECT_SOURCE_DIR}/cmake/SdpbVersion.cmake
  BYPRODUCTS "${generated_dir}/sdpb_version.hxx"
  VERBATIM)

# Settings for all SDPB targets, linked to sdpb_util
add_library(sdpb_options INTERFACE)
add_dependencies(sdpb_options sdpb_version)
target_include_directories(sdpb_options INTERFACE src external "${generated_dir}")
# -g: stack traces show source locations.
# -UNDEBUG: keep assert() enabled in dependency headers, e.g. RapidJSON, as in waf builds.
target_compile_options(sdpb_options INTERFACE -Wall -Wextra -Werror=return-type -g -UNDEBUG)
target_compile_definitions(sdpb_options INTERFACE
  ${stacktrace_define}
  # Boost.Stacktrace requires it on macOS
  $<$<PLATFORM_ID:Darwin>:_GNU_SOURCE>)

add_subdirectory(src/sdpb_util)
add_subdirectory(src/pmp)
add_subdirectory(src/pmp2sdp)
add_subdirectory(src/pmp_read)
add_subdirectory(src/sdp_solve)
add_subdirectory(src/sdpa_solve)
add_subdirectory(src/sdpb)
add_subdirectory(src/pvm2sdp)
add_subdirectory(src/sdp2input)
add_subdirectory(src/outer_limits)
add_subdirectory(src/approx_objective)
add_subdirectory(src/pmp2functions)
add_subdirectory(src/spectrum)

include(GNUInstallDirs)
install(TARGETS sdpb pmp2sdp pvm2sdp sdp2input outer_limits approx_objective pmp2functions spectrum)

enable_testing()
add_subdirectory(test)
