# Writes <OUTPUT> defining SDPB_VERSION_STRING as <VERSION>, or, if it is empty, from 'git describe'.
# Runs at every build:
#   cmake -DVERSION=<version> -DSOURCE_DIR=<dir> -DGIT_EXECUTABLE=<git> -DOUTPUT=<file> -P SdpbVersion.cmake

if(VERSION)
  set(version "${VERSION}")
else()
  set(version unknown)
  if(GIT_EXECUTABLE)
    execute_process(COMMAND "${GIT_EXECUTABLE}" describe --tags --always --dirty
      WORKING_DIRECTORY "${SOURCE_DIR}"
      RESULT_VARIABLE result OUTPUT_VARIABLE git_version OUTPUT_STRIP_TRAILING_WHITESPACE ERROR_QUIET)
    if(result EQUAL 0)
      set(version "${git_version}")
    endif()
  endif()
endif()

# Escape for a C string literal
string(REPLACE "\\" "\\\\" escaped "${version}")
string(REPLACE "\"" "\\\"" escaped "${escaped}")
set(content "#pragma once\n#define SDPB_VERSION_STRING \"${escaped}\"\n")
set(old_content "")
if(EXISTS "${OUTPUT}")
  file(READ "${OUTPUT}" old_content)
endif()
# Write only on change, so that files including it are not recompiled
if(NOT content STREQUAL old_content)
  file(WRITE "${OUTPUT}" "${content}")
  message(STATUS "SDPB version: ${version}")
endif()
