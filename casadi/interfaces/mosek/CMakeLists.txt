cmake_minimum_required(VERSION 2.8.6)
include_directories(${MOSEK_INCLUDE_DIR})

casadi_plugin(SocpSolver mosek
  mosek_socp_interface.hpp
  mosek_socp_interface.cpp
  mosek_socp_interface_meta.cpp)

casadi_plugin_link_libraries(SocpSolver mosek ${MOSEK_LIBRARIES})
