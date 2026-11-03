# Finds MPSolve and defines target MPSolve::MPSolve. Requires GMP::gmp.
# Search prefix: MPSolve_ROOT. Result variables: MPSolve_INCLUDE_DIR, MPSolve_LIBRARY, MPSolve_VERSION.

find_path(MPSolve_INCLUDE_DIR mps/mps.h)
find_library(MPSolve_LIBRARY mps)
mark_as_advanced(MPSolve_INCLUDE_DIR MPSolve_LIBRARY)

if(MPSolve_INCLUDE_DIR)
  include(SdpbHeaderVersion)
  sdpb_header_version(MPSolve_VERSION "${MPSolve_INCLUDE_DIR}/mps/version.h"
    MPS_MAJOR_VERSION MPS_MINOR_VERSION MPS_PATCH_VERSION)
endif()

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(MPSolve
  REQUIRED_VARS MPSolve_LIBRARY MPSolve_INCLUDE_DIR
  VERSION_VAR MPSolve_VERSION
  REASON_FAILURE_MESSAGE
  "Install MPSolve ${MPSolve_FIND_VERSION} or later. Set CMAKE_PREFIX_PATH or MPSolve_ROOT to its installation prefix, see Install.md.")

if(MPSolve_FOUND AND NOT TARGET MPSolve::MPSolve)
  add_library(MPSolve::MPSolve UNKNOWN IMPORTED)
  set_target_properties(MPSolve::MPSolve PROPERTIES
    IMPORTED_LOCATION "${MPSolve_LIBRARY}"
    INTERFACE_INCLUDE_DIRECTORIES "${MPSolve_INCLUDE_DIR}"
    INTERFACE_LINK_LIBRARIES GMP::gmp)
endif()
