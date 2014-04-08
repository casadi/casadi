cmake_minimum_required(VERSION 2.6)
set(CONTROL_SRCS
  dple_solver.hpp
  dple_solver.cpp
  dple_internal.hpp
  dple_internal.cpp
  simple_indef_dple_solver.hpp
  simple_indef_dple_solver.cpp
  simple_indef_dple_internal.hpp
  simple_indef_dple_internal.cpp)

casadi_library(casadi_control ${CONTROL_SRCS})

install(DIRECTORY ./
  DESTINATION include/casadi/control
  FILES_MATCHING PATTERN "*.hpp"
  PATTERN ".svn" EXCLUDE)
