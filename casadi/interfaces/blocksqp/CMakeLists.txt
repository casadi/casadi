cmake_minimum_required(VERSION 2.8.6)
include_directories(${BLOCKSQP_INCLUDE_DIR})
include_directories(${QPOASES_INCLUDE_DIR})

casadi_plugin(Nlpsol blocksqp
  blocksqp_interface.hpp
  blocksqp_interface.cpp
  blocksqp_interface_meta.cpp)

casadi_plugin_link_libraries(Nlpsol blocksqp ${BLOCKSQP_LIBRARIES} ${QPOASES_LIBRARIES})
