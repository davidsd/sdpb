# Finds MPFR and defines target MPFR::MPFR. Requires GMP::gmp.
# Search prefix: MPFR_ROOT. Result variables: MPFR_INCLUDE_DIR, MPFR_LIBRARY, MPFR_VERSION.

find_package(PkgConfig QUIET)
if(PkgConfig_FOUND)
  pkg_check_modules(PC_MPFR QUIET mpfr)
endif()
find_path(MPFR_INCLUDE_DIR mpfr.h HINTS ${PC_MPFR_INCLUDE_DIRS})
find_library(MPFR_LIBRARY mpfr HINTS ${PC_MPFR_LIBRARY_DIRS})
mark_as_advanced(MPFR_INCLUDE_DIR MPFR_LIBRARY)

if(MPFR_INCLUDE_DIR)
  include(SdpbHeaderVersion)
  sdpb_header_version(MPFR_VERSION "${MPFR_INCLUDE_DIR}/mpfr.h"
    MPFR_VERSION_MAJOR MPFR_VERSION_MINOR MPFR_VERSION_PATCHLEVEL)
endif()

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(MPFR
  REQUIRED_VARS MPFR_LIBRARY MPFR_INCLUDE_DIR
  VERSION_VAR MPFR_VERSION
  REASON_FAILURE_MESSAGE "Set CMAKE_PREFIX_PATH or MPFR_ROOT to its installation prefix, see Install.md.")

if(MPFR_FOUND AND NOT TARGET MPFR::MPFR)
  add_library(MPFR::MPFR UNKNOWN IMPORTED)
  set_target_properties(MPFR::MPFR PROPERTIES
    IMPORTED_LOCATION "${MPFR_LIBRARY}"
    INTERFACE_INCLUDE_DIRECTORIES "${MPFR_INCLUDE_DIR}"
    INTERFACE_LINK_LIBRARIES GMP::gmp)
endif()
