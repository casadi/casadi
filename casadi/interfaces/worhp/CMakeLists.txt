cmake_minimum_required(VERSION 2.6)

include_directories(../)
include_directories(${WORHP_INCLUDE_DIR})

casadi_plugin(NlpSolver worhp
  worhp_interface.hpp
  worhp_interface.cpp
  worhp_interface_meta.cpp)
casadi_plugin_link_libraries(NlpSolver worhp ${WORHP_LIBRARIES})

install(DIRECTORY ./
  DESTINATION include/casadi/interfaces/worhp
  FILES_MATCHING PATTERN "*.hpp"
  PATTERN ".svn" EXCLUDE)

if("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang" OR "${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
set_target_properties(casadi_nlpsolver_worhp PROPERTIES COMPILE_FLAGS "-Wno-unused-function")
endif()
