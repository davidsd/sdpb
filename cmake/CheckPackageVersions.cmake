# Detect packages shadowed by another installation of the same package.
#
# The build merges include directories of all dependencies in the order of SDPB_DEPENDENCIES.
# If e.g. Elemental's prefix also contains FLINT headers, they can shadow the ones found for FLINT.
# Libraries are linked by full path, but at run time the loader searches all rpath directories,
# so a library with the same soname can be shadowed as well.
# For each package, we print its header and library versions twice:
# with the package's target first ("expected"), and in the build order ("actual").

# Builds and runs a program printing "header <version>" and, if <library> is set, "library <version>".
# Sets <out> to its output, or to an empty string on failure.
function(_sdpb_version_probe out include header library)
  set(source "${include}\n#include <cstdio>\n#define STR_(x) #x\n#define STR(x) STR_(x)\nint main()\n{\n")
  string(APPEND source "  std::printf(\"header %s\\n\", ${header});\n")
  if(library)
    string(APPEND source "  std::printf(\"library %s\\n\", ${library});\n")
  endif()
  string(APPEND source "}\n")
  set(dir "${PROJECT_BINARY_DIR}/CMakeFiles/CheckPackageVersions")
  file(WRITE "${dir}/probe.cxx" "${source}")
  try_run(_sdpb_probe_run _sdpb_probe_compile "${dir}" "${dir}/probe.cxx"
    LINK_LIBRARIES ${ARGN}
    RUN_OUTPUT_VARIABLE output)
  if(_sdpb_probe_compile AND _sdpb_probe_run EQUAL 0)
    set(${out} "${output}" PARENT_SCOPE)
  else()
    set(${out} "" PARENT_SCOPE)
  endif()
endfunction()

# Files of <target> found in directories of targets preceding it in <build_order>
function(_sdpb_shadowing_files out build_order target include)
  string(REGEX MATCHALL "<[^>]+>" headers "${include}")
  string(REGEX REPLACE "[<>]" "" headers "${headers}")
  get_target_property(own_dirs ${target} INTERFACE_INCLUDE_DIRECTORIES)
  set(own_library "")
  get_target_property(type ${target} TYPE)
  if(NOT type STREQUAL "INTERFACE_LIBRARY")
    get_target_property(location ${target} LOCATION)
    get_filename_component(dir "${location}" DIRECTORY)
    list(APPEND own_dirs "${dir}")
    # e.g. libflint.so.19 -> libflint
    get_filename_component(own_library "${location}" NAME)
    string(REGEX REPLACE "\\..*" "" own_library "${own_library}")
  endif()

  list(FIND build_order ${target} index)
  list(SUBLIST build_order 0 ${index} preceding)
  set(result "")
  foreach(other IN LISTS preceding)
    get_target_property(dirs ${other} INTERFACE_INCLUDE_DIRECTORIES)
    foreach(dir IN LISTS dirs)
      if(NOT dir IN_LIST own_dirs)
        foreach(header IN LISTS headers)
          if(EXISTS "${dir}/${header}")
            list(APPEND result "${dir}/${header}")
          endif()
        endforeach()
      endif()
    endforeach()
    get_target_property(type ${other} TYPE)
    if(own_library AND NOT type STREQUAL "INTERFACE_LIBRARY")
      get_target_property(location ${other} LOCATION)
      get_filename_component(dir "${location}" DIRECTORY)
      if(NOT dir IN_LIST own_dirs)
        file(GLOB libraries "${dir}/${own_library}.*")
        list(APPEND result ${libraries})
      endif()
    endif()
  endforeach()
  list(REMOVE_DUPLICATES result)
  set(${out} "${result}" PARENT_SCOPE)
endfunction()

# Checks <target> in <build_order>. Sets <passed> to FALSE in the caller scope on failure.
function(_sdpb_check_package_version build_order target include header library)
  # e.g. FLINT::FLINT -> FLINT
  string(REGEX REPLACE "::.*" "" package ${target})
  message(CHECK_START "Checking ${package} version")
  set(others ${build_order})
  list(REMOVE_ITEM others ${target})
  _sdpb_version_probe(expected "${include}" "${header}" "${library}" ${target} ${others})
  _sdpb_version_probe(actual "${include}" "${header}" "${library}" ${build_order})
  if(NOT expected OR NOT actual)
    message(CHECK_FAIL "skipped: failed to build or run")
    set(passed FALSE PARENT_SCOPE)
    return()
  endif()
  if(expected STREQUAL actual)
    string(REGEX MATCH "header ([^\n]*)" _ "${expected}")
    message(CHECK_PASS "${CMAKE_MATCH_1}")
    return()
  endif()
  message(CHECK_FAIL "shadowed")

  set(text "The build sees a different ${package} than the one found for it:\n")
  foreach(kind IN ITEMS header library)
    string(REGEX MATCH "${kind} ([^\n]*)" _ "${expected}")
    set(expected_version "${CMAKE_MATCH_1}")
    string(REGEX MATCH "${kind} ([^\n]*)" _ "${actual}")
    if(NOT CMAKE_MATCH_1 STREQUAL expected_version)
      string(APPEND text "  ${kind}: ${CMAKE_MATCH_1}, expected ${expected_version}\n")
    endif()
  endforeach()
  string(APPEND text "Another ${package} installation shadows it:\n")
  _sdpb_shadowing_files(files "${build_order}" ${target} "${include}")
  foreach(file IN LISTS files)
    string(APPEND text "  ${file}\n")
  endforeach()
  string(APPEND text "Install ${package} in a separate directory.")
  message(SEND_ERROR "${text}")
  set(passed FALSE PARENT_SCOPE)
endfunction()

# sdpb_check_package_versions(<target or library>...)
# Checks packages in the given build order.
# After all checks pass, runs again only if include directories or library locations change.
function(sdpb_check_package_versions)
  # Plain library names, e.g. dl, have no include directories to check
  set(build_order "")
  foreach(item IN LISTS ARGN)
    if(TARGET ${item})
      list(APPEND build_order ${item})
    endif()
  endforeach()
  set(key "")
  foreach(target IN LISTS build_order)
    get_target_property(dirs ${target} INTERFACE_INCLUDE_DIRECTORIES)
    get_target_property(type ${target} TYPE)
    set(location "")
    if(NOT type STREQUAL "INTERFACE_LIBRARY")
      get_target_property(location ${target} LOCATION)
    endif()
    string(APPEND key "${target} ${dirs} ${location} | ")
  endforeach()
  if(key STREQUAL "${_SDPB_CHECKED_PACKAGE_VERSIONS}")
    return()
  endif()

  set(passed TRUE)
  _sdpb_check_package_version("${build_order}" GMP::gmpxx
    [[#include <gmp.h>]]
    [[STR(__GNU_MP_VERSION) "." STR(__GNU_MP_VERSION_MINOR) "." STR(__GNU_MP_VERSION_PATCHLEVEL)]]
    [[gmp_version]])
  _sdpb_check_package_version("${build_order}" MPFR::MPFR
    [[#include <mpfr.h>]] [[MPFR_VERSION_STRING]] [[mpfr_get_version()]])
  _sdpb_check_package_version("${build_order}" Boost::headers
    [[#include <boost/version.hpp>]] [[BOOST_LIB_VERSION]] "")
  _sdpb_check_package_version("${build_order}" Elemental::El
    [[#include <El/config.h>]] [[EL_VERSION_MAJOR "." EL_VERSION_MINOR " " EL_GIT_SHA1]] "")
  _sdpb_check_package_version("${build_order}" LibXml2::LibXml2
    [[#include <libxml/parser.h>]] [[LIBXML_DOTTED_VERSION]] [[xmlParserVersion]])
  _sdpb_check_package_version("${build_order}" RapidJSON::RapidJSON
    [[#include <rapidjson/rapidjson.h>]] [[RAPIDJSON_VERSION_STRING]] "")
  _sdpb_check_package_version("${build_order}" LibArchive::LibArchive
    [[#include <archive.h>]] [[ARCHIVE_VERSION_STRING]] [[archive_version_string()]])
  _sdpb_check_package_version("${build_order}" FLINT::FLINT
    [[#include <flint/flint.h>]] [[FLINT_VERSION]] [[flint_version]])
  # mps/version.h declares mps_get_version() without extern "C"
  _sdpb_check_package_version("${build_order}" MPSolve::MPSolve
    [[extern "C" {
#include <mps/version.h>
}]]
    [[STR(MPS_MAJOR_VERSION) "." STR(MPS_MINOR_VERSION) "." STR(MPS_PATCH_VERSION)]]
    [[mps_get_version()]])

  if(passed)
    set(_SDPB_CHECKED_PACKAGE_VERSIONS "${key}" CACHE INTERNAL "")
  endif()
endfunction()
