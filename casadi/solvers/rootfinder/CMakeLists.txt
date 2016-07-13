cmake_minimum_required(VERSION 2.8.6)

casadi_plugin(Rootfinder newton
  newton.hpp newton.cpp newton_meta.cpp)

casadi_plugin(Rootfinder nlpsol
  implicit_to_nlp.hpp implicit_to_nlp.cpp implicit_to_nlp_meta.cpp)
