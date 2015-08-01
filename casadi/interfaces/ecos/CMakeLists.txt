cmake_minimum_required(VERSION 2.8.6)
include_directories(${ECOS_INCLUDE_DIR})

casadi_plugin(SocpSolver ecos
  ecos_interface.hpp
  ecos_interface.cpp
  ecos_interface_meta.cpp)

casadi_plugin_link_libraries(SocpSolver ecos ${ECOS_LIBRARIES})
