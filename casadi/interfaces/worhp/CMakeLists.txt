cmake_minimum_required(VERSION 3.10.2)

include_directories(../)

casadi_plugin(Nlpsol worhp
  worhp_interface.hpp
  worhp_interface.cpp
  worhp_interface_meta.cpp)
casadi_plugin_link_libraries(Nlpsol worhp worhp::worhp)

