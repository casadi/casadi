cmake_minimum_required(VERSION 3.10.2)

casadi_plugin(Nlpsol uno
  uno_interface.hpp
  uno_interface.cpp
  uno_interface_meta.cpp)

casadi_plugin_link_libraries(Nlpsol uno uno)
