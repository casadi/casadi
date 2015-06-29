message(STATUS "Looking for ecos")

find_path(ECOS_INCLUDE_DIR
    ecos.h
  HINTS $ENV{ECOS}/include /usr/local/include/ecos ~/local/lib/ecos
)

if(ECOS_INCLUDE_DIR)
  message(STATUS "Found ecos include directory: ${ECOS_INCLUDE_DIR}")
else()
  message(STATUS "Could not find ecos include dir")
endif()

find_library(ECOS_LIBRARY
  NAMES ecos
  PATHS ~/local/lib /usr/local/lib
  $ENV{ECOS}
)

if(ECOS_LIBRARY)
  set(FOUND_ECOS_LIBS TRUE)
  set(ECOS_LIBRARIES ${ECOS_LIBRARY})
endif()

if(ECOS_LIBRARIES)
  set(ECOS_LIBRARIES ${ECOS_LIBRARIES})
  message(STATUS "Found ecos libraries ${ECOS_LIBRARIES}")
  set(ECOS_FOUND_LIBS TRUE)
else()
  set(ECOS_FOUND_LIBS FALSE)
  message(STATUS "Could not find ecos libraries ${ECOS_LIBRARIES}")
endif()

if(ECOS_INCLUDE_DIR AND ECOS_FOUND_LIBS)
  set(ECOS_FOUND TRUE)
else()
  set(ECOS_FOUND FALSE)
  message(STATUS "ECOS: Cound not find ecos. Try setting ECOS env var.")
endif()


