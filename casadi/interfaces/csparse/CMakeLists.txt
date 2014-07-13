cmake_minimum_required(VERSION 2.6)
include_directories(${CSPARSE_INCLUDE_DIR})

casadi_library(casadi_linearsolver_csparse
  csparse_internal.hpp
  csparse_internal.cpp)
target_link_libraries(casadi_linearsolver_csparse ${CSPARSE_LIBRARIES})

casadi_library(casadi_linearsolver_csparsecholesky
  csparse_cholesky_internal.hpp
  csparse_cholesky_internal.cpp)
target_link_libraries(casadi_linearsolver_csparsecholesky ${CSPARSE_LIBRARIES})

install(DIRECTORY ./
  DESTINATION include/casadi/interfaces/csparse
  FILES_MATCHING PATTERN "*.hpp"
  PATTERN ".svn" EXCLUDE)
