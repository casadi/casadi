cmake_minimum_required(VERSION 3.16.3)

casadi_plugin(Blas blasfeo
  blas_blasfeo.hpp
  blas_blasfeo.cpp
  blas_blasfeo_meta.cpp)
casadi_plugin_link_libraries(Blas blasfeo blasfeo)
