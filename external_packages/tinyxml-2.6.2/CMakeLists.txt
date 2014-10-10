cmake_minimum_required(VERSION 2.6)

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

if(ENABLE_STATIC)
add_library(casadi_tinyxml STATIC ${TINYXML_SRC})
endif()
if(ENABLE_SHARED)
add_library(casadi_tinyxml SHARED ${TINYXML_SRC})
endif()
install(TARGETS casadi_tinyxml
  LIBRARY DESTINATION lib
  ARCHIVE DESTINATION lib)
install(DIRECTORY ./
  DESTINATION include/casadi/external_packages/tinyxml-2.6.2
  FILES_MATCHING PATTERN "*.h"
  PATTERN ".svn" EXCLUDE)
set_target_properties(casadi_tinyxml PROPERTIES SOVERSION ${PACKAGE_SO_VERSION})
