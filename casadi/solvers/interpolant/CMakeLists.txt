cmake_minimum_required(VERSION 2.8.6)

# Linear interpolant
casadi_plugin(Interpolant linear
  linear_interpolant.hpp linear_interpolant.cpp linear_interpolant_meta.cpp
)

casadi_plugin(Interpolant bspline
  bspline_interpolant.hpp bspline_interpolant.cpp bspline_interpolant_meta.cpp
)
