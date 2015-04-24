cmake_minimum_required(VERSION 2.8.6)
include_directories(${MOSEK_INCLUDE_DIR})

casadi_plugin(SdpSolver mosek
  mosek_interface.hpp
  mosek_interface.cpp
  mosek_interface_meta.cpp)

casadi_plugin_link_libraries(SdpSolver mosek ${MOSEK_LIBRARIES})
