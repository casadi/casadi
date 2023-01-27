cmake_minimum_required(VERSION 3.10.2)

casadi_plugin(Nlpsol ampl
  ampl_interface.hpp
  ampl_interface.cpp
  ampl_interface_meta.cpp)
