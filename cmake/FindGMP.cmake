# Finds GMP with C++ support and defines targets GMP::gmp and GMP::gmpxx.
# Search prefix: GMP_ROOT. Result variables:
# GMP_INCLUDE_DIR, GMPXX_INCLUDE_DIR, GMP_LIBRARY, GMPXX_LIBRARY, GMP_VERSION.

find_package(PkgConfig QUIET)
if(PkgConfig_FOUND)
  pkg_check_modules(PC_GMPXX QUIET gmpxx)
endif()
find_path(GMP_INCLUDE_DIR gmp.h HINTS ${PC_GMPXX_INCLUDE_DIRS})
find_path(GMPXX_INCLUDE_DIR gmpxx.h HINTS ${PC_GMPXX_INCLUDE_DIRS})
find_library(GMP_LIBRARY gmp HINTS ${PC_GMPXX_LIBRARY_DIRS})
find_library(GMPXX_LIBRARY gmpxx HINTS ${PC_GMPXX_LIBRARY_DIRS})
mark_as_advanced(GMP_INCLUDE_DIR GMPXX_INCLUDE_DIR GMP_LIBRARY GMPXX_LIBRARY)

if(GMP_INCLUDE_DIR)
  include(SdpbHeaderVersion)
  sdpb_header_version(GMP_VERSION "${GMP_INCLUDE_DIR}/gmp.h"
    __GNU_MP_VERSION __GNU_MP_VERSION_MINOR __GNU_MP_VERSION_PATCHLEVEL)
endif()

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(GMP
  REQUIRED_VARS GMP_LIBRARY GMPXX_LIBRARY GMP_INCLUDE_DIR GMPXX_INCLUDE_DIR
  VERSION_VAR GMP_VERSION
  REASON_FAILURE_MESSAGE "Set CMAKE_PREFIX_PATH or GMP_ROOT to its installation prefix, see Install.md.")

if(GMP_FOUND AND NOT TARGET GMP::gmp)
  add_library(GMP::gmp UNKNOWN IMPORTED)
  set_target_properties(GMP::gmp PROPERTIES
    IMPORTED_LOCATION "${GMP_LIBRARY}"
    INTERFACE_INCLUDE_DIRECTORIES "${GMP_INCLUDE_DIR}")
  add_library(GMP::gmpxx UNKNOWN IMPORTED)
  set_target_properties(GMP::gmpxx PROPERTIES
    IMPORTED_LOCATION "${GMPXX_LIBRARY}"
    INTERFACE_INCLUDE_DIRECTORIES "${GMPXX_INCLUDE_DIR}"
    INTERFACE_LINK_LIBRARIES GMP::gmp)
endif()
