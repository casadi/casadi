cmake_minimum_required(VERSION 2.8.12)

casadi_plugin(Conic cbc
  cbc_interface.hpp
  cbc_interface.cpp
  cbc_interface_meta.cpp)
casadi_plugin_link_libraries(Conic cbc cbc)

