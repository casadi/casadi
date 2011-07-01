cmake_minimum_required(VERSION 2.0)
SET(TINYXML_SRC 
  tinystr.cpp
  tinyxml.cpp
  tinyxmlerror.cpp
  tinyxmlparser.cpp
  dump.cpp
)

if(ENABLE_STATIC)
add_library(tinyxml STATIC ${TINYXML_SRC})
endif(ENABLE_STATIC)
if(ENABLE_SHARED)
add_library(tinyxml SHARED ${TINYXML_SRC})
endif(ENABLE_SHARED)
install(TARGETS tinyxml
  LIBRARY DESTINATION lib
  ARCHIVE DESTINATION lib
)

add_executable(xmltest xmltest.cpp)
target_link_libraries(xmltest tinyxml)

add_executable(test_dump test_dump.cpp)
target_link_libraries(test_dump tinyxml)
