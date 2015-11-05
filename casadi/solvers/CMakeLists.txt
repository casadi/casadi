cmake_minimum_required(VERSION 2.8.6)

# Base classes for fixed-step integrators
casadi_library(casadi_integrators
  fixed_step_integrator.hpp
  fixed_step_integrator.cpp
  implicit_fixed_step_integrator.hpp
  implicit_fixed_step_integrator.cpp)

# Explicit Runge-Kutta integrator
casadi_plugin(Integrator rk
  rk_integrator.hpp
  rk_integrator.cpp
  rk_integrator_meta.cpp)
target_link_libraries(casadi_integrator_rk casadi_integrators)

# Collocation integrator
casadi_plugin(Integrator collocation
  collocation_integrator.hpp
  collocation_integrator.cpp
  collocation_integrator_meta.cpp)
target_link_libraries(casadi_integrator_collocation casadi_integrators)

# Old (all-at-once) collocation integrator
casadi_plugin(Integrator oldcollocation
  old_collocation_integrator.hpp
  old_collocation_integrator.cpp
  old_collocation_integrator_meta.cpp)

casadi_plugin(Rootfinder newton
  newton.hpp newton.cpp newton_meta.cpp)

# SQPMethod -  A basic SQP method
casadi_plugin(Nlpsol sqpmethod
  sqpmethod.hpp sqpmethod.cpp sqpmethod_meta.cpp)

# SCPgen -  An implementation of Lifted Newton SQP
casadi_plugin(Nlpsol scpgen
  scpgen.hpp scpgen.cpp scpgen_meta.cpp)

# Reformulations
casadi_plugin(Qpsol nlpsol
  qp_to_nlp.hpp qp_to_nlp.cpp qp_to_nlp_meta.cpp)
casadi_plugin(Rootfinder nlpsol
  implicit_to_nlp.hpp implicit_to_nlp.cpp implicit_to_nlp_meta.cpp)
casadi_plugin(LinearSolver symbolicqr
  symbolic_qr.hpp symbolic_qr.cpp symbolic_qr_meta.cpp
)

if(WITH_DL AND NOT WIN32)
  # Simple just-in-time compiler, using shell commands
  casadi_plugin(Compiler shell
    shell_compiler.hpp
    shell_compiler.cpp
    shell_compiler_meta.cpp)
endif()
