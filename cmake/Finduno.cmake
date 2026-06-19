# prefer explicit user-provided root (case-insensitive, keep both)
if(NOT DEFINED UNO_ROOT AND DEFINED uno_ROOT)
  set(UNO_ROOT "${uno_ROOT}" CACHE PATH "Root path to UNO installation" FORCE)
endif()

if(NOT DEFINED UNO_ROOT)
  set(UNO_ROOT $ENV{UNO} CACHE PATH "Root path to UNO installation (from UNO env)" )
endif()

find_path(UNO_INCLUDE_DIR
  NAMES Uno_C_API.h
  HINTS "${UNO_ROOT}/include" "${UNO_ROOT}/include/uno" "${CMAKE_PREFIX_PATH}/include"
  PATH_SUFFIXES ""
  DOC "UNO C API include directory")

find_library(UNO_LIBRARY
  NAMES uno libuno
  HINTS "${UNO_ROOT}/lib" "${UNO_ROOT}/lib64" "${CMAKE_PREFIX_PATH}/lib"
  DOC "UNO shared library")

if(UNO_INCLUDE_DIR AND UNO_LIBRARY)
  set(UNO_FOUND TRUE)
  set(UNO_INCLUDE_DIRS "${UNO_INCLUDE_DIR}")
  set(UNO_LIBRARIES "${UNO_LIBRARY}")
endif()

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(uno DEFAULT_MSG UNO_LIBRARIES UNO_INCLUDE_DIRS)
mark_as_advanced(UNO_ROOT UNO_INCLUDE_DIR UNO_LIBRARY)