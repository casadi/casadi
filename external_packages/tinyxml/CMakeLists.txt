cmake_minimum_required(VERSION 2.6)
set(TINYXML_SRC
  tinystr.cpp
  tinyxml.cpp
  tinyxmlerror.cpp
  tinyxmlparser.cpp
  dump.cpp
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
  DESTINATION include/casadi/external_packages/tinyxml
  FILES_MATCHING PATTERN "*.h"
  PATTERN ".svn" EXCLUDE)

add_executable(test_dump test_dump.cpp)
target_link_libraries(test_dump ${TINYXML_LIBRARIES})
