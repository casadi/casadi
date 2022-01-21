cmake_minimum_required(VERSION 2.8.6)

casadi_plugin(Conic clp
  clp_interface.hpp
  clp_interface.cpp
  clp_interface_meta.cpp)
casadi_plugin_link_libraries(Conic clp clp)
