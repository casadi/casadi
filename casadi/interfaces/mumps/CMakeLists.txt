cmake_minimum_required(VERSION 2.8.6)
include_directories(${MUMPS_INCLUDE_DIR})

casadi_plugin(Linsol mumps
  mumps_interface.hpp
  mumps_interface.cpp
  mumps_interface_meta.cpp)
casadi_plugin_link_libraries(Linsol mumps ${MUMPS_LIBRARIES})
