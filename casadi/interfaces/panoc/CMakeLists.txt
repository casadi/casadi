cmake_minimum_required(VERSION 2.8.6)
include_directories(${PANOC_INCLUDE_DIR})

message("Include dir: ${PANOC_INCLUDE_DIR}")
set(NLPSOL_PANOC_SRCS
  panoc.cpp
  panoc.hpp
  panoc_meta.cpp)

casadi_plugin(Nlpsol panoc ${NLPSOL_PANOC_SRCS})
casadi_plugin_link_libraries(Nlpsol panoc ${PANOC_LIBRARIES})
