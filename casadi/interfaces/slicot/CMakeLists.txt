cmake_minimum_required(VERSION 2.6)

casadi_library(casadi_slicot_interface
  slicot_tools.hpp
  slicot_tools.cpp)
target_link_libraries(casadi_slicot_interface casadi_control ${SLICOT_LIBRARIES})
set_property(GLOBAL APPEND PROPERTY CASADI_SWIG_MODULES casadi_slicot_interface)

casadi_library(casadi_dplesolver_slicot
  psd_indef_dple_internal.hpp
  psd_indef_dple_internal.cpp)
target_link_libraries(casadi_dplesolver_slicot casadi_slicot_interface casadi_control ${SLICOT_LIBRARIES})

install(DIRECTORY ./
  DESTINATION include/casadi/interfaces/slicot
  FILES_MATCHING PATTERN "*.hpp"
  PATTERN ".svn" EXCLUDE)

