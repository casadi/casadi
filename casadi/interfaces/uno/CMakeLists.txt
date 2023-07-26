cmake_minimum_required(VERSION 3.10.2)

casadi_plugin(Nlpsol uno
  uno_interface.hpp
  uno_interface.cpp
  uno_interface_meta.cpp)

find_library(ma57 NAMES ma57 HINTS $ENV{MA57})
# target_link_libraries(uno INTERFACE ${ma57})
casadi_plugin_link_libraries(Nlpsol uno uno ma57)
