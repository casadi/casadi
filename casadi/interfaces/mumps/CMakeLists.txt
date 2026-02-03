cmake_minimum_required(VERSION 3.16.3)

casadi_plugin(Linsol mumps
  mumps_interface.hpp
  mumps_interface.cpp
  mumps_interface_meta.cpp)
casadi_plugin_link_libraries(Linsol mumps mumps)
