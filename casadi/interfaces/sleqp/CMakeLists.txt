cmake_minimum_required(VERSION 3.10.2)

set(NLPSOL_SLEQP_SRCS
  sleqp_interface.cpp
  sleqp_func.cpp
)

if(WITH_SLEQP_CALLBACK)
  add_definitions(-DWITH_SLEQP_CALLBACK)
endif()

casadi_plugin(Nlpsol sleqp ${NLPSOL_SLEQP_SRCS})
casadi_plugin_link_libraries(Nlpsol sleqp sleqp::sleqp)
