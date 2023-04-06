cmake_minimum_required(VERSION 3.10.2)

casadi_plugin(Nlpsol snopt
  snopt_interface.hpp
  snopt_interface.cpp
  snopt_interface_meta.cpp)

casadi_plugin_link_libraries(Nlpsol snopt snopt::snopt)
