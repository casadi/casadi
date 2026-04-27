cmake_minimum_required(VERSION 3.10.2)

casadi_plugin(Blas classic
  blas_classic.hpp
  blas_classic.cpp
  blas_classic_meta.cpp)
casadi_plugin_link_libraries(Blas classic lapack)
