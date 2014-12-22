cmake_minimum_required(VERSION 2.8.6)
include_directories(${DSDP_INCLUDE_DIR})

casadi_plugin(SdpSolver sdp dsdp
  dsdp_interface.hpp
  dsdp_interface.cpp
  dsdp_interface_meta.cpp)

casadi_plugin_link_libraries(sdp dsdp ${DSDP_LIBRARIES})
