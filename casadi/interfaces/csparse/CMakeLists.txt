cmake_minimum_required(VERSION 2.6)
include_directories(${CSPARSE_INCLUDE_DIR})

set(CSPARSE_INTERFACE_SRCS
  csparse.hpp
  csparse.cpp
  csparse_internal.hpp
  csparse_internal.cpp
  csparse_cholesky.hpp
  csparse_cholesky.cpp
  csparse_cholesky_internal.hpp
  csparse_cholesky_internal.cpp)

casadi_library(casadi_csparse_interface ${CSPARSE_INTERFACE_SRCS})

target_link_libraries(casadi_csparse_interface ${CSPARSE_LIBRARIES})

install(DIRECTORY ./
  DESTINATION include/casadi/interfaces/csparse
  FILES_MATCHING PATTERN "*.hpp"
  PATTERN ".svn" EXCLUDE)
