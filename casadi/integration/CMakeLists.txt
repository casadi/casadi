cmake_minimum_required(VERSION 2.6)
casadi_library_swig(casadi_integration
  fixed_step_integrator_internal.hpp
  fixed_step_integrator_internal.cpp
  implicit_fixed_step_integrator_internal.hpp
  implicit_fixed_step_integrator_internal.cpp
  integration_tools.cpp
  integration_tools.hpp
)

# Explicit Runge-Kutta integrator
casadi_library(casadi_integrator_rk
  rk_integrator_internal.hpp
  rk_integrator_internal.cpp)
target_link_libraries(casadi_integrator_rk casadi_integration)

# Collocation integrator
casadi_library(casadi_integrator_collocation
  collocation_integrator_internal.hpp
  collocation_integrator_internal.cpp)
target_link_libraries(casadi_integrator_collocation casadi_integration)

# Old (all-at-once) collocation integrator
casadi_library(casadi_integrator_oldcollocation
  old_collocation_integrator_internal.hpp
  old_collocation_integrator_internal.cpp)
target_link_libraries(casadi_integrator_oldcollocation casadi_integration)

install(DIRECTORY ./
  DESTINATION include/casadi/integration
  FILES_MATCHING PATTERN "*.hpp"
  PATTERN ".svn" EXCLUDE)
