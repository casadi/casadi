cmake_minimum_required(VERSION 2.6)
include_directories(${CSPARSE_INCLUDE_DIR})

casadi_plugin(LinearSolver csparse
  csparse_interface.hpp
  csparse_interface.cpp
  csparse_interface_meta.cpp)
casadi_plugin_link_libraries(LinearSolver csparse ${CSPARSE_LIBRARIES})

casadi_plugin(LinearSolver csparsecholesky
  csparse_cholesky_internal.hpp
  csparse_cholesky_internal.cpp
  csparse_cholesky_internal_meta.cpp
  )
casadi_plugin_link_libraries(LinearSolver csparsecholesky ${CSPARSE_LIBRARIES})

install(DIRECTORY ./
  DESTINATION include/casadi/interfaces/csparse
  FILES_MATCHING PATTERN "*.hpp"
  PATTERN ".svn" EXCLUDE)
