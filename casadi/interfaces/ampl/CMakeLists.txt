cmake_minimum_required(VERSION 3.16.3)

casadi_plugin(Nlpsol ampl
  ampl_interface.hpp
  ampl_interface.cpp
  ampl_interface_meta.cpp)
