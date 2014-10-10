cmake_minimum_required(VERSION 2.6)
include_directories(${QPOASES_INCLUDE_DIR})

casadi_plugin(QpSolver qpoases
  qpoases_interface.hpp
  qpoases_interface.cpp
  qpoases_interface_meta.cpp)

casadi_plugin_link_libraries(QpSolver qpoases ${QPOASES_LIBRARIES})

install(DIRECTORY ./
  DESTINATION include/casadi/interfaces/qpoases
  FILES_MATCHING PATTERN "*.hpp"
  PATTERN ".svn" EXCLUDE)

if("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang" OR "${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
set_target_properties(casadi_qpsolver_qpoases PROPERTIES COMPILE_FLAGS "-Wno-overloaded-virtual")
endif()
