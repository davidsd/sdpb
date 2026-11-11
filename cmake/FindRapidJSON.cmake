# Finds the header-only RapidJSON library and defines target RapidJSON::RapidJSON.
# Search prefix: RapidJSON_ROOT. Result variables: RapidJSON_INCLUDE_DIR, RapidJSON_VERSION.

find_package(PkgConfig QUIET)
if(PkgConfig_FOUND)
  pkg_check_modules(PC_RapidJSON QUIET RapidJSON)
endif()
find_path(RapidJSON_INCLUDE_DIR rapidjson/rapidjson.h HINTS ${PC_RapidJSON_INCLUDE_DIRS})
mark_as_advanced(RapidJSON_INCLUDE_DIR)

if(RapidJSON_INCLUDE_DIR)
  include(SdpbHeaderVersion)
  sdpb_header_version(RapidJSON_VERSION "${RapidJSON_INCLUDE_DIR}/rapidjson/rapidjson.h"
    RAPIDJSON_MAJOR_VERSION RAPIDJSON_MINOR_VERSION RAPIDJSON_PATCH_VERSION)
endif()

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(RapidJSON
  REQUIRED_VARS RapidJSON_INCLUDE_DIR
  VERSION_VAR RapidJSON_VERSION
  REASON_FAILURE_MESSAGE "Set CMAKE_PREFIX_PATH or RapidJSON_ROOT to its installation prefix, see Install.md.")

if(RapidJSON_FOUND AND NOT TARGET RapidJSON::RapidJSON)
  add_library(RapidJSON::RapidJSON INTERFACE IMPORTED)
  set_target_properties(RapidJSON::RapidJSON PROPERTIES
    INTERFACE_INCLUDE_DIRECTORIES "${RapidJSON_INCLUDE_DIR}")
endif()
