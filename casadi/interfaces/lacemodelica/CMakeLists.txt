cmake_minimum_required(VERSION 3.10.2)

casadi_plugin(ModelicaParser lacemodelica
  lace_modelica.hpp
  lace_modelica.cpp
  lace_modelica_meta.cpp)
casadi_plugin_link_libraries(ModelicaParser lacemodelica lacemodelica::lacemodelica)
