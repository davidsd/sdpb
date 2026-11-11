# Finds a BLAS library with the C interface (cblas.h), OpenBLAS by default,
# and defines target CBLAS::CBLAS.
# Search prefix: CBLAS_ROOT. For another BLAS, set CBLAS_LIBRARY to its full path.
# Result variables: CBLAS_INCLUDE_DIR, CBLAS_LIBRARY.

find_package(PkgConfig QUIET)
if(PkgConfig_FOUND)
  pkg_check_modules(PC_CBLAS QUIET openblas)
endif()
find_path(CBLAS_INCLUDE_DIR cblas.h HINTS ${PC_CBLAS_INCLUDE_DIRS} PATH_SUFFIXES openblas)
find_library(CBLAS_LIBRARY openblas HINTS ${PC_CBLAS_LIBRARY_DIRS})
mark_as_advanced(CBLAS_INCLUDE_DIR CBLAS_LIBRARY)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(CBLAS
  REQUIRED_VARS CBLAS_LIBRARY CBLAS_INCLUDE_DIR
  REASON_FAILURE_MESSAGE
  "Set CMAKE_PREFIX_PATH or CBLAS_ROOT to the OpenBLAS installation prefix, or CBLAS_LIBRARY to another BLAS library. See Install.md.")

if(CBLAS_FOUND AND NOT TARGET CBLAS::CBLAS)
  add_library(CBLAS::CBLAS UNKNOWN IMPORTED)
  set_target_properties(CBLAS::CBLAS PROPERTIES
    IMPORTED_LOCATION "${CBLAS_LIBRARY}"
    INTERFACE_INCLUDE_DIRECTORIES "${CBLAS_INCLUDE_DIR}")
endif()
