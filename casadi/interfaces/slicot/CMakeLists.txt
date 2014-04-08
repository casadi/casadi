cmake_minimum_required(VERSION 2.6)

set(SLICOT_INTERFACE_SRCS
  slicot_tools.hpp
  slicot_tools.cpp
  psd_indef_dple_solver.hpp
  psd_indef_dple_solver.cpp
  psd_indef_dple_internal.hpp
  psd_indef_dple_internal.cpp)

casadi_library(casadi_slicot_interface ${SLICOT_INTERFACE_SRCS})
target_link_libraries(casadi_slicot_interface casadi_control ${SLICOT_LIBRARIES})

install(DIRECTORY ./
  DESTINATION include/casadi/interfaces/slicot
  FILES_MATCHING PATTERN "*.hpp"
  PATTERN ".svn" EXCLUDE)

