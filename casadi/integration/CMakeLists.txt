cmake_minimum_required(VERSION 2.6)
set(INTEGRATION_SRCS
  fixed_step_integrator.hpp
  fixed_step_integrator.cpp
  fixed_step_integrator_internal.hpp
  fixed_step_integrator_internal.cpp
  implicit_fixed_step_integrator.hpp
  implicit_fixed_step_integrator.cpp
  implicit_fixed_step_integrator_internal.hpp
  implicit_fixed_step_integrator_internal.cpp
  rk_integrator.hpp
  rk_integrator.cpp
  rk_integrator_internal.hpp
  rk_integrator_internal.cpp
  collocation_integrator.hpp
  collocation_integrator.cpp
  collocation_integrator_internal.hpp
  collocation_integrator_internal.cpp
  old_collocation_integrator.hpp
  old_collocation_integrator.cpp
  old_collocation_integrator_internal.hpp
  old_collocation_integrator_internal.cpp
  integration_tools.cpp
  integration_tools.hpp
)

casadi_library(casadi_integration ${INTEGRATION_SRCS})

install(DIRECTORY ./
  DESTINATION include/casadi/integration
  FILES_MATCHING PATTERN "*.hpp"
  PATTERN ".svn" EXCLUDE)
