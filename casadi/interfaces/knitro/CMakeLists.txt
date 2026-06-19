cmake_minimum_required(VERSION 3.16.3)

casadi_plugin(Nlpsol knitro
  knitro_interface.hpp
  knitro_interface.cpp
  knitro_interface_meta.cpp)

casadi_plugin_link_libraries(Nlpsol knitro knitro::knitro)
