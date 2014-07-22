cmake_minimum_required(VERSION 2.6)
include_directories(${CPLEX_INCLUDE_DIR})

casadi_plugin(QpSolver cplex
  cplex_internal.hpp
  cplex_internal.cpp
  cplex_internal_meta.cpp)
casadi_plugin_link_libraries(QpSolver cplex ${CPLEX_LIBRARIES})

install(DIRECTORY ./
  DESTINATION include/casadi/interfaces/cplex
  FILES_MATCHING PATTERN "*.hpp"
  PATTERN ".svn" EXCLUDE)


