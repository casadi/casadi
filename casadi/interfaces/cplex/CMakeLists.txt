cmake_minimum_required(VERSION 3.10.2)

casadi_plugin(Conic cplex
  cplex_interface.hpp
  cplex_interface.cpp
  cplex_interface_meta.cpp)
casadi_plugin_link_libraries(Conic cplex cplex::cplex)

if(CPLEX_DEFINITIONS)
  add_definitions(${CPLEX_DEFINITIONS})
endif()
