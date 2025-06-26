message(STATUS "Looking for MUMPS")

find_path(MUMPS_INCLUDE_DIR
    mumps_mpi.h "mumps_seq/mpi.h"
    HINTS $ENV{MUMPS}/include
    PATH_SUFFIXES "coin-or/mumps"
)

if(MUMPS_INCLUDE_DIR)
  message(STATUS "Found MUMPS include directory: ${MUMPS_INCLUDE_DIR}")
else()
  message(STATUS "Could not find MUMPS include dir")
endif()

# libraries
set(MUMPS_LIBS_LIST
  dmumps_seq coinmumps dmumps
)

set(MUMPS_LIBRARIES)
foreach(LIB ${MUMPS_LIBS_LIST})
  find_library(MUMPS_LIB_${LIB}
    NAMES ${LIB}
    HINTS $ENV{MUMPS}/lib)
  if(MUMPS_LIB_${LIB})
#    message(STATUS "Found ${LIB}: ${MUMPS_LIB_${LIB}}")
    set(MUMPS_LIBRARIES ${MUMPS_LIB_${LIB}})
#  else()
#    message(STATUS "Could not find lib${LIB}")
    break()
  endif()
endforeach()

if(MUMPS_LIBRARIES)
  set(MUMPS_LIBRARIES ${MUMPS_LIBRARIES})
  message(STATUS "Found MUMPS libraries ${MUMPS_LIBRARIES}")
  set(MUMPS_FOUND_LIBS TRUE)
else()
  message(STATUS "Could not find MUMPS libraries")
endif()

if(MUMPS_INCLUDE_DIR AND MUMPS_FOUND_LIBS)
  set(MUMPS_FOUND TRUE)
else()
  message(STATUS "MUMPS: Cound not find MUMPS. Try setting MUMPS env var.")
endif()

if(MUMPS_FOUND)
  add_library(mumps INTERFACE)
  target_link_libraries(mumps INTERFACE ${MUMPS_LIBRARIES})
  target_include_directories(mumps INTERFACE ${MUMPS_INCLUDE_DIRS})
endif()
