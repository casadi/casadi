cmake_minimum_required(VERSION 2.0)
SET(TINYXML_SRC 
  tinystr.cpp
  tinyxml.cpp
  tinyxmlerror.cpp
  tinyxmlparser.cpp
  dump.cpp
)

add_library(tinyxml STATIC ${TINYXML_SRC})
if(ENABLE_SHARED)
add_library(tinyxml_shared SHARED ${TINYXML_SRC})
endif(ENABLE_SHARED)

add_executable(xmltest xmltest.cpp)
target_link_libraries(xmltest tinyxml)

add_executable(test_dump test_dump.cpp)
target_link_libraries(test_dump tinyxml)
