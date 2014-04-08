cmake_minimum_required(VERSION 2.6)

include_directories(../)
include_directories(${WORHP_INCLUDE_DIR})

set(WORHP_INTERFACE_SRCS
  worhp_internal.hpp
  worhp_internal.cpp
  worhp_solver.hpp
  worhp_solver.cpp
)

casadi_library(casadi_worhp_interface ${WORHP_INTERFACE_SRCS})
target_link_libraries(casadi_worhp_interface ${WORHP_LIBRARIES})

install(DIRECTORY ./
  DESTINATION include/casadi/interfaces/worhp
  FILES_MATCHING PATTERN "*.hpp"
  PATTERN ".svn" EXCLUDE)

if("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang" OR "${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
set_target_properties(casadi_worhp_interface PROPERTIES COMPILE_FLAGS "-Wno-unused-function")
endif()
