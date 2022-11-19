cmake_minimum_required(VERSION 2.8.6)


casadi_plugin(Conic highs
  highs_interface.hpp
  highs_interface.cpp
  highs_interface_meta.cpp)
casadi_plugin_link_libraries(Conic highs highs)
