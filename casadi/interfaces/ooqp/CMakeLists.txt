cmake_minimum_required(VERSION 3.10.2)

include_directories(${OOQP_INCLUDE_DIRS})

casadi_plugin(Conic ooqp
  ooqp_interface.hpp
  ooqp_interface.cpp
  ooqp_interface_meta.cpp)

casadi_plugin_link_libraries(Conic ooqp ${OOQP_LIBRARIES})
