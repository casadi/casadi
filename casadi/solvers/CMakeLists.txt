cmake_minimum_required(VERSION 2.6)

# Explicit Runge-Kutta integrator
casadi_library(casadi_integrator_rk
  fixed_step_integrator_internal.hpp
  fixed_step_integrator_internal.cpp
  rk_integrator_internal.hpp
  rk_integrator_internal.cpp)

# Collocation integrator
casadi_library(casadi_integrator_collocation
  implicit_fixed_step_integrator_internal.hpp
  implicit_fixed_step_integrator_internal.cpp
  collocation_integrator_internal.hpp
  collocation_integrator_internal.cpp)
target_link_libraries(casadi_integrator_collocation casadi_integrator_rk)

# Old (all-at-once) collocation integrator
casadi_library(casadi_integrator_oldcollocation
  old_collocation_integrator_internal.hpp
  old_collocation_integrator_internal.cpp)

casadi_library(casadi_implicitfunction_newton
  newton_implicit_internal.hpp   newton_implicit_internal.cpp)

casadi_library(casadi_homotopynlpsolver_simple
   simple_homotopy_nlp_internal.hpp   simple_homotopy_nlp_internal.cpp)

# SQPMethod -  A basic SQP method
casadi_library(casadi_nlpsolver_sqpmethod   sqp_internal.hpp   sqp_internal.cpp)

# SCPgen -  An implementation of Lifted Newton SQP
casadi_library(casadi_nlpsolver_scpgen   scpgen_internal.hpp   scpgen_internal.cpp)

# StabilizedSQP -  An implementation of Stabilized SQP
casadi_library(casadi_nlpsolver_stabilizedsqp   stabilized_sqp_internal.hpp   stabilized_sqp_internal.cpp)

# Reformulations
casadi_library(casadi_qpsolver_nlp nlp_qp_internal.hpp   nlp_qp_internal.cpp)
casadi_library(casadi_implicitfunction_nlp   nlp_implicit_internal.hpp   nlp_implicit_internal.cpp)
casadi_library(casadi_lpsolver_qp   qp_lp_internal.cpp   qp_lp_internal.hpp)
casadi_library(casadi_qpsolver_qcqp   qcqp_qp_internal.cpp   qcqp_qp_internal.hpp)
casadi_library(casadi_socpsolver_sdp   sdp_socp_internal.cpp   sdp_socp_internal.hpp)
casadi_library(casadi_stabilizedqpsolver_qp   qp_stabilizer_internal.cpp   qp_stabilizer_internal.hpp)
casadi_library(casadi_linearsolver_symbolicqr   symbolic_qr_internal.hpp   symbolic_qr_internal.cpp)
if(WITH_CSPARSE)
  casadi_library(casadi_qcqpsolver_socp   socp_qcqp_internal.cpp   socp_qcqp_internal.hpp)
  target_link_libraries(casadi_qcqpsolver_socp   casadi_linearsolver_csparsecholesky)
  casadi_library(casadi_sdqpsolver_sdp   sdp_sdqp_internal.cpp    sdp_sdqp_internal.hpp)
  target_link_libraries(casadi_sdqpsolver_sdp   casadi_linearsolver_csparsecholesky)
endif()

install(DIRECTORY ./
  DESTINATION include/casadi/solvers
  FILES_MATCHING PATTERN "*.hpp"
  PATTERN ".svn" EXCLUDE)
