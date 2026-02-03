cmake_minimum_required(VERSION 3.16.3)

casadi_plugin(Nlpsol snopt
  snopt_interface.hpp
  snopt_interface.cpp
  snopt_interface_meta.cpp)

casadi_plugin_link_libraries(Nlpsol snopt snopt::snopt)
