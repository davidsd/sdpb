# Finds Elemental (https://gitlab.com/bootstrapcollaboration/elemental)
# and defines target Elemental::El. Requires MPI::MPI_CXX, GMP::gmpxx and MPFR::MPFR.
# Search prefix: Elemental_ROOT. Result variables: Elemental_INCLUDE_DIR, Elemental_LIBRARY,
# Elemental_pmrrr_LIBRARY, Elemental_ElSuiteSparse_LIBRARY, Elemental_VERSION.
#
# We do not use ElementalConfig.cmake installed by Elemental:
# its target El lacks Elemental's include directory,
# and find_package() can pick a stale Elemental build tree from the CMake package registry.

find_path(Elemental_INCLUDE_DIR El.hpp)
find_library(Elemental_LIBRARY El)
find_library(Elemental_pmrrr_LIBRARY pmrrr)
find_library(Elemental_ElSuiteSparse_LIBRARY ElSuiteSparse)
mark_as_advanced(Elemental_INCLUDE_DIR Elemental_LIBRARY Elemental_pmrrr_LIBRARY Elemental_ElSuiteSparse_LIBRARY)

# Version is e.g. "0.88-dev"
if(Elemental_INCLUDE_DIR AND EXISTS "${Elemental_INCLUDE_DIR}/El/config.h")
  file(STRINGS "${Elemental_INCLUDE_DIR}/El/config.h" lines REGEX "^#define EL_VERSION_(MAJOR|MINOR) ")
  if(lines MATCHES "EL_VERSION_MAJOR \"([^\"]*)\".*EL_VERSION_MINOR \"([^\"]*)\"")
    set(Elemental_VERSION "${CMAKE_MATCH_1}.${CMAKE_MATCH_2}")
  endif()
endif()

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(Elemental
  REQUIRED_VARS Elemental_LIBRARY Elemental_pmrrr_LIBRARY Elemental_ElSuiteSparse_LIBRARY Elemental_INCLUDE_DIR
  VERSION_VAR Elemental_VERSION
  REASON_FAILURE_MESSAGE "Set CMAKE_PREFIX_PATH or Elemental_ROOT to its installation prefix, see Install.md.")

if(Elemental_FOUND AND NOT TARGET Elemental::El)
  add_library(Elemental::El UNKNOWN IMPORTED)
  set_target_properties(Elemental::El PROPERTIES
    IMPORTED_LOCATION "${Elemental_LIBRARY}"
    INTERFACE_INCLUDE_DIRECTORIES "${Elemental_INCLUDE_DIR}"
    INTERFACE_LINK_LIBRARIES
    "${Elemental_pmrrr_LIBRARY};${Elemental_ElSuiteSparse_LIBRARY};MPI::MPI_CXX;GMP::gmpxx;MPFR::MPFR")
endif()
