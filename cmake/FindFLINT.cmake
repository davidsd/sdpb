# Finds FLINT and defines target FLINT::FLINT. Requires GMP::gmp and MPFR::MPFR.
# Search prefix: FLINT_ROOT. Result variables: FLINT_INCLUDE_DIR, FLINT_LIBRARY, FLINT_VERSION.

find_package(PkgConfig QUIET)
if(PkgConfig_FOUND)
  pkg_check_modules(PC_FLINT QUIET flint)
endif()
find_path(FLINT_INCLUDE_DIR flint/flint.h HINTS ${PC_FLINT_INCLUDE_DIRS})
find_library(FLINT_LIBRARY flint HINTS ${PC_FLINT_LIBRARY_DIRS})
mark_as_advanced(FLINT_INCLUDE_DIR FLINT_LIBRARY)

if(FLINT_INCLUDE_DIR)
  include(SdpbHeaderVersion)
  sdpb_header_version(FLINT_VERSION "${FLINT_INCLUDE_DIR}/flint/flint.h"
    __FLINT_VERSION __FLINT_VERSION_MINOR __FLINT_VERSION_PATCHLEVEL)
endif()

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(FLINT
  REQUIRED_VARS FLINT_LIBRARY FLINT_INCLUDE_DIR
  VERSION_VAR FLINT_VERSION
  REASON_FAILURE_MESSAGE
  "Install FLINT ${FLINT_FIND_VERSION} or later. Set CMAKE_PREFIX_PATH or FLINT_ROOT to its installation prefix, see Install.md.")

if(FLINT_FOUND AND NOT TARGET FLINT::FLINT)
  add_library(FLINT::FLINT UNKNOWN IMPORTED)
  set_target_properties(FLINT::FLINT PROPERTIES
    IMPORTED_LOCATION "${FLINT_LIBRARY}"
    INTERFACE_INCLUDE_DIRECTORIES "${FLINT_INCLUDE_DIR}"
    INTERFACE_LINK_LIBRARIES "GMP::gmp;MPFR::MPFR")
endif()
