cmake_minimum_required(VERSION 3.10.2)

casadi_plugin(Linsol ma27
  ma27_interface.hpp
  ma27_interface.cpp
  ma27_interface_meta.cpp)
casadi_plugin_link_libraries(Linsol ma27 hsl::hsl)
