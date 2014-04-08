# TRY TO FIND THE INCLUDE DIRECTORY
find_library(SLICOT_LIB
  slicot
  PATHS $ENV{SLICOT_LIBRARY_DIR})

if(NOT SLICOT_LIB)
  message(STATUS "SLICOT: Cound not find library slicot. Try stetting SLICOT_LIBRARY_DIR env var.")
endif()

if(SLICOT_LIB)
  set(SLICOT_LIBRARIES ${SLICOT_LIB})
  set(SLICOT_FOUND TRUE)
endif()
