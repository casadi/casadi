cmake_minimum_required(VERSION 2.8.6)

# Explicit Runge-Kutta integrator
casadi_plugin(Integrator rk
  rk_integrator.hpp
  rk_integrator.cpp
  rk_integrator_meta.cpp)

# Collocation integrator
casadi_plugin(Integrator collocation
  collocation_integrator.hpp
  collocation_integrator.cpp
  collocation_integrator_meta.cpp)
