cmake_minimum_required(VERSION 2.6)
include_directories(${DSDP_INCLUDE_DIR})

set(DSDP_INTERFACE_SRCS
  dsdp_internal.hpp
  dsdp_internal.cpp
  dsdp_solver.hpp
  dsdp_solver.cpp)

casadi_library(casadi_dsdp_interface ${DSDP_INTERFACE_SRCS})
target_link_libraries(casadi_dsdp_interface ${DSDP_LIBRARIES})

install(DIRECTORY ./
  DESTINATION include/casadi/interfaces/dsdp
  FILES_MATCHING PATTERN "*.hpp"
  PATTERN ".svn" EXCLUDE)
