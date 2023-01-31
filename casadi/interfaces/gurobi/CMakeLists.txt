cmake_minimum_required(VERSION 3.10.2)

include_directories(../)
include_directories(${GUROBI_INCLUDE_DIRS})

casadi_plugin(Conic gurobi
  gurobi_interface.hpp
  gurobi_interface.cpp
  gurobi_interface_meta.cpp)

casadi_plugin_link_libraries(Conic gurobi ${GUROBI_LIBRARIES})
