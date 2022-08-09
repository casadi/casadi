cmake_minimum_required(VERSION 2.8.6)

casadi_plugin(Conic hpipm
  hpipm_interface.hpp
  hpipm_interface.cpp
  hpipm_interface_meta.cpp)

casadi_plugin_link_libraries(Conic hpipm hpipm)