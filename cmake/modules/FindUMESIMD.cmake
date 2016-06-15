# - Locate UMESIMD libraries and includes
# Defines:
#
#  UMESIMD_FOUND
#  UMESIMD_INCLUDE_DIR

if(NOT UMESIMD_DIR)
   find_path(UMESIMD_DIR "umesimd/UMESimd.h")
endif()

find_path(UMESIMD_INCLUDE_DIR umesimd/UMESimd.h
          HINTS ${UMESIMD_DIR})

message(STATUS "UMESimd includes path: ${UMESIMD_INCLUDE_DIR}")

INCLUDE(FindPackageHandleStandardArgs)
FIND_PACKAGE_HANDLE_STANDARD_ARGS(UMESIMD DEFAULT_MSG UMESIMD_INCLUDE_DIR)
mark_as_advanced(UMESIMD_FOUND UMESIMD_INCLUDE_DIR)
