cmake_minimum_required(VERSION 2.8.6)

include_directories(../)

casadi_plugin(Conic superscs
  superscs_interface.hpp
  superscs_interface.cpp
  superscs_interface_meta.cpp)

casadi_plugin_link_libraries(Conic superscs superscs)
