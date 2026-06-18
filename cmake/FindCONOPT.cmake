# Prefer explicit user-provided root (case-insensitive, keep both)
if(NOT DEFINED CONOPT_ROOT AND DEFINED conopt_ROOT)
  set(CONOPT_ROOT "${conopt_ROOT}" CACHE PATH "Root path to CONOPT installation" FORCE)
endif()

if(NOT DEFINED CONOPT_ROOT)
  set(CONOPT_ROOT $ENV{CONOPT} CACHE PATH "Root path to CONOPT installation (from CONOPT env)" )
endif()

# Find the include directory containing conopt.h
find_path(CONOPT_INCLUDE_DIR
  NAMES conopt.h
  HINTS "${CONOPT_ROOT}/include" "${CONOPT_ROOT}" "${CMAKE_PREFIX_PATH}/include"
  PATH_SUFFIXES ""
  DOC "CONOPT C API include directory")

if(CONOPT_INCLUDE_DIR)
   message(STATUS "Found CONOPT include dir: ${CONOPT_INCLUDE_DIR}")
else()
   message(STATUS "Could not find CONOPT include dir")
endif()

# Find the CONOPT library
find_library(CONOPT_LIBRARY
  NAMES conopt libconopt
  HINTS "${CONOPT_ROOT}/lib" "${CONOPT_ROOT}/lib64" "${CONOPT_ROOT}" "${CMAKE_PREFIX_PATH}/lib"
  DOC "CONOPT library")

if(CONOPT_LIBRARY)
  set(CONOPT_LIBRARIES "${CONOPT_LIBRARY}")
  message(STATUS "Found CONOPT libraries: ${CONOPT_LIBRARIES}")
else()
  message(STATUS "Could not find CONOPT library")
endif()

# Handle standard arguments and set CONOPT_FOUND
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(CONOPT DEFAULT_MSG CONOPT_LIBRARIES CONOPT_INCLUDE_DIR)

# Create the imported target so casadi_plugin_link_libraries(Nlpsol conopt conopt::conopt) works
if(CONOPT_FOUND)
  if(NOT TARGET conopt::conopt)
    add_library(conopt::conopt UNKNOWN IMPORTED)
    set_target_properties(conopt::conopt PROPERTIES
      INTERFACE_INCLUDE_DIRECTORIES "${CONOPT_INCLUDE_DIR}"
      IMPORTED_LOCATION "${CONOPT_LIBRARIES}"
    )
  endif()
endif()

# Hide these variables from the standard CMake GUI to keep it clean
mark_as_advanced(CONOPT_ROOT CONOPT_INCLUDE_DIR CONOPT_LIBRARY)
