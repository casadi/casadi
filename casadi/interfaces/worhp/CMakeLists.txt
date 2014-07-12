cmake_minimum_required(VERSION 2.6)

include_directories(../)
include_directories(${WORHP_INCLUDE_DIR})

casadi_library(casadi_nlpsolver_worhp
  worhp_internal.hpp
  worhp_internal.cpp)
target_link_libraries(casadi_nlpsolver_worhp ${WORHP_LIBRARIES})

install(DIRECTORY ./
  DESTINATION include/casadi/interfaces/worhp
  FILES_MATCHING PATTERN "*.hpp"
  PATTERN ".svn" EXCLUDE)

if("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang" OR "${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
set_target_properties(casadi_nlpsolver_worhp PROPERTIES COMPILE_FLAGS "-Wno-unused-function")
endif()
