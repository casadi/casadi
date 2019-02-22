cmake_minimum_required(VERSION 2.8.6)

include_directories(${CBC_INCLUDE_DIRS})
link_directories(${CBC_LIBRARY_DIRS})

casadi_plugin(Conic cbc
  cbc_interface.hpp
  cbc_interface.cpp
  cbc_interface_meta.cpp)
casadi_plugin_link_libraries(Conic cbc ${CBC_LIBRARIES})
