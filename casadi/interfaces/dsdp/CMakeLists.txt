cmake_minimum_required(VERSION 2.6)
include_directories(${DSDP_INCLUDE_DIR})

casadi_plugin(SdpSolver dsdp
  dsdp_interface.hpp
  dsdp_interface.cpp
  dsdp_interface_meta.cpp)

casadi_plugin_link_libraries(SdpSolver dsdp ${DSDP_LIBRARIES})

install(DIRECTORY ./
  DESTINATION include/casadi/interfaces/dsdp
  FILES_MATCHING PATTERN "*.hpp"
  PATTERN ".svn" EXCLUDE)
