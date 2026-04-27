cmake_minimum_required(VERSION 3.10.2)

casadi_plugin(Blas blasfeo
  blas_blasfeo.hpp
  blas_blasfeo.cpp
  blas_blasfeo_meta.cpp)
casadi_plugin_link_libraries(Blas blasfeo blasfeo)
