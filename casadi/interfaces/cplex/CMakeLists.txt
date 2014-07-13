cmake_minimum_required(VERSION 2.6)
include_directories(${CPLEX_INCLUDE_DIR})

casadi_library(casadi_qpsolver_cplex
  cplex_internal.hpp
  cplex_internal.cpp)

target_link_libraries(casadi_qpsolver_cplex ${CPLEX_LIBRARIES})

install(DIRECTORY ./
  DESTINATION include/casadi/interfaces/cplex
  FILES_MATCHING PATTERN "*.hpp"
  PATTERN ".svn" EXCLUDE)


