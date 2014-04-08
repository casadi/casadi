cmake_minimum_required(VERSION 2.6)
include_directories(${CPLEX_INCLUDE_DIR})

set(CPLEX_INTERFACE_SRCS
  cplex_solver.hpp
  cplex_solver.cpp
  cplex_internal.hpp
  cplex_internal.cpp)

casadi_library(casadi_cplex_interface ${CPLEX_INTERFACE_SRCS})
target_link_libraries(casadi_cplex_interface ${CPLEX_LIBRARIES})

install(DIRECTORY ./
  DESTINATION include/casadi/interfaces/cplex
  FILES_MATCHING PATTERN "*.hpp"
  PATTERN ".svn" EXCLUDE)


