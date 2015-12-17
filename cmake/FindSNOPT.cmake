# C headers
find_path(SNOPT_INCLUDE_DIR
snopt_cwrap.h
HINTS $ENV{SNOPT}/include
)

if(SNOPT_INCLUDE_DIR)
   set(SNOPT_FOUND_INCLUDE TRUE)
   message(STATUS "Found SNOPT include dir: ${SNOPT_INCLUDE_DIR}")
else()
   message(STATUS "Could not find SNOPT include dir")
endif()

# Library
find_library(SNOPT_LIBRARY_C
snopt7
HINTS $ENV{SNOPT}/lib
)

# Library
find_library(SNOPT_LIBRARY
snopt7_c
HINTS $ENV{SNOPT}/lib
)

if(SNOPT_LIBRARY AND SNOPT_LIBRARY_C)
  set(SNOPT_LIBRARIES ${SNOPT_LIBRARY} ${SNOPT_LIBRARY_C})
  message(STATUS "Found SNOPT library: ${SNOPT_LIBRARY}")
else()
  message(STATUS "Could not find SNOPT libs")
endif()

if(SNOPT_FOUND_INCLUDE AND SNOPT_LIBRARIES)
  set(SNOPT_FOUND TRUE)
endif()
