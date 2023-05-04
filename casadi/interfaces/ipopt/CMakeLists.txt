cmake_minimum_required(VERSION 3.10.2)

set(NLPSOL_IPOPT_SRCS
  ipopt_interface.hpp
  ipopt_interface.cpp
  ipopt_nlp.hpp
  ipopt_nlp.cpp
  ipopt_interface_meta.cpp)

if(WITH_IPOPT_CALLBACK)
  add_definitions(-DWITH_IPOPT_CALLBACK)
endif()

casadi_plugin(Nlpsol ipopt ${NLPSOL_IPOPT_SRCS})
casadi_plugin_link_libraries(Nlpsol ipopt ipopt)

#set_property(TARGET casadi_nlpsol_ipopt PROPERTY PKG_CONFIG_REQUIRES ipopt)
