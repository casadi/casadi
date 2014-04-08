cmake_minimum_required(VERSION 2.6)
include_directories(${KNITRO_INCLUDE_DIR})

set(KNITRO_INTERFACE_SRCS
  knitro_solver.hpp
  knitro_solver.cpp
  knitro_internal.hpp
  knitro_internal.cpp)

casadi_library(casadi_knitro_interface ${KNITRO_INTERFACE_SRCS})
target_link_libraries(casadi_knitro_interface ${KNITRO_LIBRARIES})

install(DIRECTORY ./
  DESTINATION include/casadi/interfaces/knitro
  FILES_MATCHING PATTERN "*.hpp"
  PATTERN ".svn" EXCLUDE)
