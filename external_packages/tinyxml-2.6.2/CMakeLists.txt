cmake_minimum_required(VERSION 2.8.6)

# cmake_policy use MACOSX_RPATH=ON by default
if(COMMAND cmake_policy)
  if(CMAKE_MAJOR_VERSION GREATER 2)
    cmake_policy(SET CMP0042 NEW)
  endif()
endif()

set(TINYXML_SRC
  tinystr.cpp
  tinyxml.cpp
  tinyxmlerror.cpp
  tinyxmlparser.cpp
)

casadi_external_library(casadi_tinyxml ${TINYXML_SRC})

install(DIRECTORY ./
  DESTINATION include/casadi/external_packages/tinyxml-2.6.2
  FILES_MATCHING PATTERN "*.h"
  PATTERN ".svn" EXCLUDE)
