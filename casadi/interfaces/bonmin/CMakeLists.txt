cmake_minimum_required(VERSION 2.8.12)
include_directories(SYSTEM ${BONMIN_INCLUDE_DIRS})
link_directories(${BONMIN_LIBRARY_DIRS})

set(NLPSOL_BONMIN_SRCS
  bonmin_interface.hpp
  bonmin_interface.cpp
  bonmin_nlp.hpp
  bonmin_nlp.cpp
  bonmin_interface_meta.cpp)


casadi_plugin(Nlpsol bonmin ${NLPSOL_BONMIN_SRCS})
casadi_plugin_link_libraries(Nlpsol bonmin bonmin)

#set_property(TARGET casadi_nlpsol_bonmin PROPERTY PKG_CONFIG_REQUIRES bonmin)
