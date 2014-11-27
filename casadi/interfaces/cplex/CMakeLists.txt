cmake_minimum_required(VERSION 2.8.6)
include_directories(${CPLEX_INCLUDE_DIR})

casadi_plugin(QpSolver cplex
  cplex_interface.hpp
  cplex_interface.cpp
  cplex_interface_meta.cpp)
casadi_plugin_link_libraries(QpSolver cplex ${CPLEX_LIBRARIES})
