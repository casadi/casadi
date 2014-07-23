cmake_minimum_required(VERSION 2.6)
casadi_library(casadi_control
  dple_solver.hpp
  dple_solver.cpp
  dple_internal.hpp
  dple_internal.cpp)
set_property(GLOBAL APPEND PROPERTY CASADI_SWIG_MODULES casadi_control)

casadi_plugin(DpleSolver simple
  simple_indef_dple_internal.hpp
  simple_indef_dple_internal.cpp
  simple_indef_dple_internal_meta.cpp)
target_link_libraries(casadi_dplesolver_simple casadi_control)

install(DIRECTORY ./
  DESTINATION include/casadi/control
  FILES_MATCHING PATTERN "*.hpp"
  PATTERN ".svn" EXCLUDE)
