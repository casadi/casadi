cmake_minimum_required(VERSION 3.16.3)

casadi_plugin(Nlpsol conopt
  conopt_interface.hpp
  conopt_interface.cpp
)

casadi_plugin_link_libraries(Nlpsol conopt conopt::conopt)
