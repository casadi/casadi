cmake_minimum_required(VERSION 2.8.6)

# Explicit Runge-Kutta integrator
casadi_plugin(Ivpsol rk
  rk_integrator.hpp
  rk_integrator.cpp
  rk_integrator_meta.cpp)

# Collocation integrator
casadi_plugin(Ivpsol collocation
  collocation_integrator.hpp
  collocation_integrator.cpp
  collocation_integrator_meta.cpp)

casadi_plugin(Nlsol newton
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
casadi_plugin(Nlsol nlpsol
  implicit_to_nlp.hpp implicit_to_nlp.cpp implicit_to_nlp_meta.cpp)
casadi_plugin(Linsol symbolicqr
  symbolic_qr.hpp symbolic_qr.cpp symbolic_qr_meta.cpp
)

if(WITH_DL AND NOT WIN32)
  # Simple just-in-time compiler, using shell commands
  casadi_plugin(Compiler shell
    shell_compiler.hpp
    shell_compiler.cpp
    shell_compiler_meta.cpp)
endif()
