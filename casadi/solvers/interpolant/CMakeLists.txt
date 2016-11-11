cmake_minimum_required(VERSION 2.8.6)

# Linear interpolant
casadi_plugin(Interpolant linear
  linear_interpolant.hpp linear_interpolant.cpp linear_interpolant_meta.cpp
)
