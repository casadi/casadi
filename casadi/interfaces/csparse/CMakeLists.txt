cmake_minimum_required(VERSION 2.8.6)
include_directories(${CSPARSE_INCLUDE_DIR})

casadi_plugin(LinearSolver linsol csparse
  csparse_interface.hpp
  csparse_interface.cpp
  csparse_interface_meta.cpp)
casadi_plugin_link_libraries(linsol csparse ${CSPARSE_LIBRARIES})

casadi_plugin(LinearSolver linsol csparsecholesky
  csparse_cholesky_internal.hpp
  csparse_cholesky_internal.cpp
  csparse_cholesky_internal_meta.cpp
  )
casadi_plugin_link_libraries(linsol csparsecholesky ${CSPARSE_LIBRARIES})
